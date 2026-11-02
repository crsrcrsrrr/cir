# Build cir~ for Pure Data
PD_INCLUDE ?= /usr/include/pd
CFLAGS ?= -O2 -Wall -Wextra -Wno-unused-parameter
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  EXT = pd_darwin
  LDFLAGS += -bundle -undefined dynamic_lookup
else
  EXT = pd_linux
  LDFLAGS += -shared
endif

cir~.$(EXT): cir~.c
	$(CC) $(CFLAGS) -fPIC -I$(PD_INCLUDE) $(LDFLAGS) -o $@ $< -lm

# Windows 64-bit, cross-compiled with zig. pd.def lists the Pd functions
# we use, so no pd.dll is needed to build.
windows: cir~.dll
cir~.dll: cir~.c pd.def
	zig dlltool -d pd.def -l libpd_imp.a -m i386:x86-64
	zig cc -target x86_64-windows-gnu -O2 -DMSW -DPD_LONGINTTYPE="long long" \
		-I$(PD_INCLUDE) -shared -s -o $@ $< libpd_imp.a
	rm -f libpd_imp.a cir~.lib cir~.pdb

clean:
	rm -f cir~.pd_linux cir~.pd_darwin cir~.dll
.PHONY: clean windows
