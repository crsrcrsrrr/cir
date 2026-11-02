/*
 * cir~ : a norns softcut-style circular buffer voice for Pure Data.
 *
 * Reads from / writes to a Pd array. Several cir~ objects pointing at the
 * same array behave like several softcut voices sharing one buffer.
 *
 * Each voice has two sub-heads. A loop wrap or a "position" jump starts the
 * idle head at the new spot and crossfades between the two over fade_time,
 * so cuts are click-free. Recording uses the same fades:
 *     buf = buf * pre_level + in * rec_level
 * weighted by each head's fade, and the write is resampled so any rate
 * (including negative and fractional) fills the buffer without gaps.
 *
 * Times and positions are in seconds, like softcut.
 */

#include "m_pd.h"
#include <math.h>
#include <string.h>

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifdef _WIN32
#define CIR_EXPORT __declspec(dllexport)
#else
#define CIR_EXPORT __attribute__((visibility("default")))
#endif

#define MAX_SPLAT 64 /* most buffer samples one head may write per tick */

static t_class *cir_tilde_class;

enum { HEAD_STOPPED, HEAD_FADEIN, HEAD_PLAY, HEAD_FADEOUT };

typedef struct _head {
    double phase;     /* position in samples */
    double prevphase; /* position on the previous tick (for writing) */
    float fade;       /* 0..1 */
    int state;
} t_head;

typedef struct _cir_tilde {
    t_object x_obj;
    t_float x_f;

    t_symbol *x_arrayname;
    t_word *x_vec;
    int x_npoints;
    t_clock *x_clock;
    int x_clockset;

    double x_sr;
    t_head x_head[2];
    int x_active;

    int x_play;
    int x_rec;
    int x_loop;
    double x_loopstart; /* seconds */
    double x_loopend;   /* seconds, < 0 = end of buffer */
    double x_fadetime;  /* seconds */
    float x_fadeinc;

    /* slewed parameters: current, target, slew time, coefficient */
    float x_rate, x_rate_t, x_rate_coef;
    double x_rate_slew;
    float x_level, x_level_t, x_level_coef;
    double x_level_slew;
    float x_reclev, x_reclev_t, x_prelev, x_prelev_t, x_recpre_coef;
    double x_recpre_slew;
    float x_recgate; /* ramps 0..1 over fade_time when rec toggles */

    float x_prev_in;
} t_cir_tilde;

/* ---------- helpers ---------- */

static inline long wrapi(long i, long n)
{
    i %= n;
    return i < 0 ? i + n : i;
}

static inline double wrapd(double p, double n)
{
    p = fmod(p, n);
    return p < 0 ? p + n : p;
}

static float slew_coef(double t, double sr)
{
    return (t > 0 && sr > 0) ? (float)exp(-1.0 / (t * sr)) : 0.f;
}

static void cir_update_coefs(t_cir_tilde *x)
{
    double ft = x->x_fadetime * x->x_sr;
    x->x_fadeinc = ft > 1 ? (float)(1.0 / ft) : 1.f;
    x->x_rate_coef = slew_coef(x->x_rate_slew, x->x_sr);
    x->x_level_coef = slew_coef(x->x_level_slew, x->x_sr);
    x->x_recpre_coef = slew_coef(x->x_recpre_slew, x->x_sr);
}

static int cir_getarray(t_cir_tilde *x, int complain)
{
    t_garray *a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class);
    x->x_vec = 0;
    x->x_npoints = 0;
    if (!a) {
        if (complain && *x->x_arrayname->s_name)
            pd_error(x, "cir~: %s: no such array", x->x_arrayname->s_name);
        return 0;
    }
    if (!garray_getfloatwords(a, &x->x_npoints, &x->x_vec)) {
        pd_error(x, "cir~: %s: bad template", x->x_arrayname->s_name);
        x->x_vec = 0;
        x->x_npoints = 0;
        return 0;
    }
    garray_usedindsp(a);
    return 1;
}

static void cir_redraw_tick(t_cir_tilde *x)
{
    t_garray *a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class);
    x->x_clockset = 0;
    if (a) garray_redraw(a);
}

/* loop bounds in samples, sanitised against the buffer size */
static void cir_loop_bounds(t_cir_tilde *x, double *start, double *end)
{
    double n = x->x_npoints;
    double s = x->x_loopstart * x->x_sr;
    double e = x->x_loopend < 0 ? n : x->x_loopend * x->x_sr;
    if (s < 0) s = 0;
    if (s > n - 2) s = n - 2;
    if (e > n) e = n;
    if (e < s + 2) e = s + 2;
    *start = s;
    *end = e;
}

/* jump: fade out the active head, fade in the other one at pos (samples) */
static void cir_cut(t_cir_tilde *x, double pos)
{
    t_head *old = &x->x_head[x->x_active];
    t_head *nu = &x->x_head[!x->x_active];
    if (old->state != HEAD_STOPPED) old->state = HEAD_FADEOUT;
    nu->phase = nu->prevphase = pos;
    nu->fade = 0;
    nu->state = HEAD_FADEIN;
    x->x_active = !x->x_active;
}

static inline float hermite(float f, float y0, float y1, float y2, float y3)
{
    float c1 = 0.5f * (y2 - y0);
    float c2 = y0 - 2.5f * y1 + 2.f * y2 - 0.5f * y3;
    float c3 = 0.5f * (y3 - y0) + 1.5f * (y1 - y2);
    return ((c3 * f + c2) * f + c1) * f + y1;
}

static inline float cir_peek(t_word *vec, long n, double phase)
{
    double fl = floor(phase);
    long i = (long)fl;
    float f = (float)(phase - fl);
    return hermite(f,
        vec[wrapi(i - 1, n)].w_float, vec[wrapi(i, n)].w_float,
        vec[wrapi(i + 1, n)].w_float, vec[wrapi(i + 2, n)].w_float);
}

/* Write the stretch of buffer the head crossed since the last tick,
   interpolating the input across it. */
static inline int cir_poke(t_word *vec, long n, t_head *h,
    float in0, float in1, float rec, float pre)
{
    double p0 = h->prevphase, p1 = h->phase, d = p1 - p0;
    long k, kend;
    int count = 0;
    if (d == 0) return 0;
    if (fabs(d) > MAX_SPLAT) {
        long i = wrapi((long)floor(p1 + 0.5), n);
        vec[i].w_float = vec[i].w_float * pre + in1 * rec;
        return 1;
    }
    if (d > 0) { k = (long)floor(p0) + 1; kend = (long)floor(p1); }
    else       { k = (long)ceil(p1);      kend = (long)ceil(p0) - 1; }
    for (; k <= kend; k++, count++) {
        float t = (float)((k - p0) / d);
        float v = in0 + (in1 - in0) * t;
        long i = wrapi(k, n);
        vec[i].w_float = vec[i].w_float * pre + v * rec;
    }
    return count;
}

/* ---------- dsp ---------- */

static t_int *cir_tilde_perform(t_int *w)
{
    t_cir_tilde *x = (t_cir_tilde *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    t_sample *out = (t_sample *)(w[3]);
    t_sample *phout = (t_sample *)(w[4]);
    int nblock = (int)(w[5]);
    t_word *vec = x->x_vec;
    long n = x->x_npoints;
    double start, end, len;
    int i, h, wrote = 0;
    float gateinc = x->x_fadeinc;

    if (!vec || n < 4) {
        memset(out, 0, sizeof(t_sample) * nblock);
        memset(phout, 0, sizeof(t_sample) * nblock);
        return (w + 6);
    }
    cir_loop_bounds(x, &start, &end);
    len = end - start;

    for (i = 0; i < nblock; i++) {
        float xin = in[i]; /* read first: in and out may share memory */
        float y = 0;
        t_head *act;

        /* parameter slews */
        x->x_rate = x->x_rate_t + (x->x_rate - x->x_rate_t) * x->x_rate_coef;
        x->x_level = x->x_level_t + (x->x_level - x->x_level_t) * x->x_level_coef;
        x->x_reclev = x->x_reclev_t + (x->x_reclev - x->x_reclev_t) * x->x_recpre_coef;
        x->x_prelev = x->x_prelev_t + (x->x_prelev - x->x_prelev_t) * x->x_recpre_coef;
        if (x->x_rec) { if ((x->x_recgate += gateinc) > 1) x->x_recgate = 1; }
        else if ((x->x_recgate -= gateinc) < 0) x->x_recgate = 0;

        for (h = 0; h < 2; h++) {
            t_head *hd = &x->x_head[h];
            if (hd->state == HEAD_STOPPED) continue;

            /* read (before write, so we hear the old material) */
            y += cir_peek(vec, n, hd->phase) * sinf(hd->fade * (float)M_PI_2);

            /* write */
            if (x->x_recgate > 0) {
                float g = x->x_recgate * hd->fade;
                wrote |= cir_poke(vec, n, hd, x->x_prev_in, xin,
                    x->x_reclev * g, 1.f - g * (1.f - x->x_prelev));
            }

            /* advance */
            hd->prevphase = hd->phase;
            hd->phase += x->x_rate;
            if (hd->state == HEAD_FADEIN) {
                if ((hd->fade += x->x_fadeinc) >= 1) { hd->fade = 1; hd->state = HEAD_PLAY; }
            } else if (hd->state == HEAD_FADEOUT) {
                if ((hd->fade -= x->x_fadeinc) <= 0) { hd->fade = 0; hd->state = HEAD_STOPPED; }
            }
        }

        /* loop / end-of-region handling on the active head */
        act = &x->x_head[x->x_active];
        if (act->state == HEAD_FADEIN || act->state == HEAD_PLAY) {
            if (x->x_rate > 0 && act->phase >= end) {
                if (x->x_loop) cir_cut(x, start + fmod(act->phase - end, len));
                else act->state = HEAD_FADEOUT;
            } else if (x->x_rate < 0 && act->phase < start) {
                if (x->x_loop) cir_cut(x, end - fmod(start - act->phase, len));
                else act->state = HEAD_FADEOUT;
            }
        }

        out[i] = y * x->x_level;
        phout[i] = (t_sample)(wrapd(x->x_head[x->x_active].phase, n) / x->x_sr);
        x->x_prev_in = xin;
    }

    if (wrote && !x->x_clockset) {
        x->x_clockset = 1;
        clock_delay(x->x_clock, 200);
    }
    return (w + 6);
}

static void cir_tilde_dsp(t_cir_tilde *x, t_signal **sp)
{
    if (sp[0]->s_sr != x->x_sr) {
        x->x_sr = sp[0]->s_sr;
        cir_update_coefs(x);
    }
    cir_getarray(x, 1);
    dsp_add(cir_tilde_perform, 5, x,
        sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, (t_int)sp[0]->s_n);
}

/* ---------- messages ---------- */

static void cir_set(t_cir_tilde *x, t_symbol *s)
{
    x->x_arrayname = s;
    cir_getarray(x, 1);
}

static void cir_play(t_cir_tilde *x, t_floatarg f)
{
    int on = f != 0;
    t_head *act = &x->x_head[x->x_active];
    if (on && (act->state == HEAD_STOPPED || act->state == HEAD_FADEOUT))
        cir_cut(x, act->phase);
    else if (!on && act->state != HEAD_STOPPED)
        act->state = HEAD_FADEOUT;
    x->x_play = on;
}

static void cir_position(t_cir_tilde *x, t_floatarg f)
{
    double p = f * x->x_sr;
    if (x->x_npoints > 0) p = wrapd(p, x->x_npoints);
    if (x->x_play) cir_cut(x, p);
    else {
        /* not playing: just move the parked head */
        t_head *act = &x->x_head[x->x_active];
        act->phase = act->prevphase = p;
    }
}

static void cir_rec(t_cir_tilde *x, t_floatarg f) { x->x_rec = f != 0; }
static void cir_loop(t_cir_tilde *x, t_floatarg f) { x->x_loop = f != 0; }
static void cir_loop_start(t_cir_tilde *x, t_floatarg f) { x->x_loopstart = f < 0 ? 0 : f; }
static void cir_loop_end(t_cir_tilde *x, t_floatarg f) { x->x_loopend = f; }

static void cir_rate(t_cir_tilde *x, t_floatarg f) { x->x_rate_t = f; }
static void cir_level(t_cir_tilde *x, t_floatarg f) { x->x_level_t = f; }
static void cir_rec_level(t_cir_tilde *x, t_floatarg f) { x->x_reclev_t = f; }
static void cir_pre_level(t_cir_tilde *x, t_floatarg f) { x->x_prelev_t = f; }

static void cir_fade_time(t_cir_tilde *x, t_floatarg f)
{
    x->x_fadetime = f < 0 ? 0 : f;
    cir_update_coefs(x);
}
static void cir_rate_slew_time(t_cir_tilde *x, t_floatarg f)
{
    x->x_rate_slew = f < 0 ? 0 : f;
    cir_update_coefs(x);
}
static void cir_level_slew_time(t_cir_tilde *x, t_floatarg f)
{
    x->x_level_slew = f < 0 ? 0 : f;
    cir_update_coefs(x);
}
static void cir_recpre_slew_time(t_cir_tilde *x, t_floatarg f)
{
    x->x_recpre_slew = f < 0 ? 0 : f;
    cir_update_coefs(x);
}

static void cir_clear_range(t_cir_tilde *x, long from, long count)
{
    long i;
    if (!x->x_vec) cir_getarray(x, 1);
    if (!x->x_vec) return;
    for (i = 0; i < count; i++)
        x->x_vec[wrapi(from + i, x->x_npoints)].w_float = 0;
    cir_redraw_tick(x);
}

static void cir_buffer_clear(t_cir_tilde *x)
{
    if (!x->x_vec) cir_getarray(x, 1);
    if (x->x_vec) cir_clear_range(x, 0, x->x_npoints);
}

static void cir_buffer_clear_region(t_cir_tilde *x, t_floatarg start, t_floatarg dur)
{
    if (!x->x_vec) cir_getarray(x, 1);
    if (!x->x_vec || dur <= 0) return;
    long count = (long)(dur * x->x_sr);
    if (count > x->x_npoints) count = x->x_npoints;
    cir_clear_range(x, (long)(start * x->x_sr), count);
}

static void cir_print(t_cir_tilde *x)
{
    post("cir~ [%s] %d samples @ %g Hz", x->x_arrayname->s_name, x->x_npoints, x->x_sr);
    post("  play %d  rec %d  loop %d  loop %g..%g s  fade_time %g s",
        x->x_play, x->x_rec, x->x_loop, x->x_loopstart,
        x->x_loopend < 0 ? (x->x_sr > 0 ? x->x_npoints / x->x_sr : 0) : x->x_loopend,
        x->x_fadetime);
    post("  rate %g  level %g  rec_level %g  pre_level %g",
        x->x_rate_t, x->x_level_t, x->x_reclev_t, x->x_prelev_t);
    post("  position %g s", x->x_sr > 0 ? x->x_head[x->x_active].phase / x->x_sr : 0);
}

/* ---------- lifecycle ---------- */

static void *cir_tilde_new(t_symbol *s)
{
    t_cir_tilde *x = (t_cir_tilde *)pd_new(cir_tilde_class);
    outlet_new(&x->x_obj, &s_signal); /* audio */
    outlet_new(&x->x_obj, &s_signal); /* position in seconds */
    x->x_arrayname = s;
    x->x_clock = clock_new(x, (t_method)cir_redraw_tick);
    x->x_sr = sys_getsr();
    if (x->x_sr <= 0) x->x_sr = 44100;

    x->x_loop = 1;
    x->x_loopstart = 0;
    x->x_loopend = -1;
    x->x_fadetime = 0.01;
    x->x_rate = x->x_rate_t = 1;
    x->x_level = x->x_level_t = 1;
    x->x_reclev = x->x_reclev_t = 1;
    x->x_prelev = x->x_prelev_t = 0;
    x->x_head[0].state = x->x_head[1].state = HEAD_STOPPED;
    cir_update_coefs(x);
    return x;
}

static void cir_tilde_free(t_cir_tilde *x)
{
    clock_free(x->x_clock);
}

CIR_EXPORT void cir_tilde_setup(void)
{
    cir_tilde_class = class_new(gensym("cir~"),
        (t_newmethod)cir_tilde_new, (t_method)cir_tilde_free,
        sizeof(t_cir_tilde), CLASS_DEFAULT, A_DEFSYM, 0);
    CLASS_MAINSIGNALIN(cir_tilde_class, t_cir_tilde, x_f);
    class_addmethod(cir_tilde_class, (t_method)cir_tilde_dsp, gensym("dsp"), A_CANT, 0);

    class_addmethod(cir_tilde_class, (t_method)cir_set, gensym("set"), A_SYMBOL, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_play, gensym("play"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_rec, gensym("rec"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_loop, gensym("loop"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_loop_start, gensym("loop_start"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_loop_end, gensym("loop_end"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_position, gensym("position"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_rate, gensym("rate"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_level, gensym("level"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_rec_level, gensym("rec_level"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_pre_level, gensym("pre_level"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_fade_time, gensym("fade_time"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_rate_slew_time, gensym("rate_slew_time"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_level_slew_time, gensym("level_slew_time"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_recpre_slew_time, gensym("recpre_slew_time"), A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_buffer_clear, gensym("buffer_clear"), 0);
    class_addmethod(cir_tilde_class, (t_method)cir_buffer_clear_region, gensym("buffer_clear_region"), A_FLOAT, A_FLOAT, 0);
    class_addmethod(cir_tilde_class, (t_method)cir_print, gensym("print"), 0);
}
