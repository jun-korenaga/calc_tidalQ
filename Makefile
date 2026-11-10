CXXFLAGS ?= -O3 -Wall -Wextra -Wpedantic
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
PKG_CONFIG ?= pkg-config

ifeq ($(shell uname -s),Darwin)
LAPACK_PKG ?= openblas
openblas_prefix := $(shell brew --prefix openblas 2>/dev/null)
ifneq ($(openblas_prefix),)
lapack_pkg_config := PKG_CONFIG_PATH="$(openblas_prefix)/lib/pkgconfig:$$PKG_CONFIG_PATH" $(PKG_CONFIG)
endif
else
LAPACK_PKG ?= lapacke
endif
lapack_pkg_config ?= $(PKG_CONFIG)
LAPACK_CPPFLAGS ?= $(patsubst -I%,-isystem %,$(shell $(lapack_pkg_config) --cflags-only-I $(LAPACK_PKG)))
LAPACK_LIBS ?= $(shell $(lapack_pkg_config) --libs $(LAPACK_PKG))

program := calc_tidalQ
sources := calc_tidalQ.cc tidal_sys_relaxS.cc earth_propS.cc util.cc error.cc
objects := $(sources:%.cc=build/%.o)

.PHONY: all install clean
all: $(program)

$(program): $(objects)
	$(CXX) $(LDFLAGS) -o $@ $(objects) $(LAPACK_LIBS) $(LDLIBS)

build/tidal_sys_relaxS.o: override CPPFLAGS += $(LAPACK_CPPFLAGS)

build/%.o: %.cc | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 -MMD -MP -c $< -o $@

build:
	mkdir -p $@

install: $(program)
	$(INSTALL) -d "$(DESTDIR)$(BINDIR)"
	$(INSTALL) -m 755 $(program) "$(DESTDIR)$(BINDIR)/$(program)"

clean:
	$(RM) -r build $(program)

-include $(objects:.o=.d)
