# Makefile: build and a versioned, relocatable install of ltr-forage and its library driver.
# After rmblast's (../rmblast-port/rmblast/Makefile).
#
# Layout produced by `make install`:
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/ltr-forage          the compiled binary
#   $(PREFIX)/$(NAME)-$(VERSION)/bin/ltr_library.py      link to ../pipeline/ltr_library.py
#   $(PREFIX)/$(NAME)-$(VERSION)/pipeline/*              the library driver and its step scripts
#   $(PREFIX)/$(NAME)-$(VERSION)/data/*                  protein profiles, scoring matrices
#   $(PREFIX)/$(NAME)-$(VERSION)/examples/*              Arabidopsis chromosome 4, the README's example
#   $(PREFIX)/$(NAME)-$(VERSION)/{README.md,CHANGELOG.md,LICENSE,docs/pipeline_figure.svg,docs/evidence_example.png,docs/accuracy-history.md}
#
# ltr-forage finds data/TEdomains.hmm relative to itself (../data beside bin/), and the driver finds its step
# scripts and data/ relative to pipeline/, so the tree can be moved and nothing needs editing after unpacking.
#
# Examples:
#   make                          # build (release)
#   make install PREFIX=/opt      # -> /opt/ltr-forage-<version>/...
#   make dist                     # a versioned tarball (staged under ./stage)
#   make version                  # print the version from Cargo.toml

PREFIX   ?= /usr/local
DESTDIR  ?=
NAME     ?= ltr-forage
CARGO    ?= cargo
INSTALL  ?= install
BIN       = ltr-forage
MANIFEST  = Cargo.toml
# Shipped files are listed by name or pattern, so that a working note dropped into the tree does not end up
# in a release by accident.
DOCS      = README.md CHANGELOG.md LICENSE
PIPELINE  = $(wildcard pipeline/*.py pipeline/*.sh) pipeline/evidence_template.html
DATA      = data/TEdomains.hmm data/TEfam.hmm
DOMAINS   = data/domains/README.md data/domains/profiles.tsv data/domains/build_domains.py
MATRICES  = $(wildcard data/matrices/*)
EXAMPLES  = examples/ath_chr4.fa.gz examples/README.md

# The version of record: [package] version in Cargo.toml.
VERSION := $(shell sed -n '/^\[package\]/,/^\[/s/^version *= *"\([^"]*\)".*/\1/p' $(MANIFEST) | head -1)
$(if $(VERSION),,$(error could not read version from $(MANIFEST)))

PKGDIR = $(DESTDIR)$(PREFIX)/$(NAME)-$(VERSION)

.PHONY: all build test version check-release install uninstall dist clean

all: build

version:
	@echo $(VERSION)

build:
	$(CARGO) build --release --bin $(BIN)

test:
	$(CARGO) test --release

# Before tagging a release (see the script for what it checks):
#   make check-release TAG=X.Y.Z
check-release: build
	bash .github/scripts/check_release.sh $(VERSION) "$(TAG)" $(BIN)

install: build
	$(INSTALL) -d $(PKGDIR)/bin $(PKGDIR)/pipeline $(PKGDIR)/data/domains $(PKGDIR)/data/matrices $(PKGDIR)/examples $(PKGDIR)/docs
	$(INSTALL) -m 0755 target/release/$(BIN) $(PKGDIR)/bin/$(BIN)
	$(INSTALL) -m 0755 $(filter %.py %.sh,$(PIPELINE)) $(PKGDIR)/pipeline/
	$(INSTALL) -m 0644 pipeline/evidence_template.html $(PKGDIR)/pipeline/
	ln -sf ../pipeline/ltr_library.py $(PKGDIR)/bin/ltr_library.py
	$(INSTALL) -m 0644 $(DATA) $(PKGDIR)/data/
	$(INSTALL) -m 0644 $(DOMAINS) $(PKGDIR)/data/domains/
	$(INSTALL) -m 0644 $(MATRICES) $(PKGDIR)/data/matrices/
	$(INSTALL) -m 0644 $(EXAMPLES) $(PKGDIR)/examples/
	$(INSTALL) -m 0644 docs/pipeline_figure.svg docs/evidence_example.png docs/accuracy-history.md $(PKGDIR)/docs/
	$(INSTALL) -m 0644 $(DOCS) $(PKGDIR)/
	@echo "Installed $(NAME) $(VERSION) -> $(PKGDIR)"
	@echo "  $$($(PKGDIR)/bin/$(BIN) --version)"

uninstall:
	rm -rf $(PKGDIR)

# Versioned tarball <name>-<version>.tar.gz of the installed hierarchy, built from a staging install: under
# DESTDIR if given (the release workflow passes one), else ./stage.  It never installs into PREFIX itself.
DIST_DESTDIR = $(if $(DESTDIR),$(DESTDIR),$(CURDIR)/stage)
dist: build
	$(MAKE) install DESTDIR=$(DIST_DESTDIR)
	tar -C $(DIST_DESTDIR)$(PREFIX) -czf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	@echo "Wrote $(NAME)-$(VERSION).tar.gz (staged in $(DIST_DESTDIR)$(PREFIX))"

clean:
	$(CARGO) clean
