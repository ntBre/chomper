run:
	cargo run

testflags =

ifdef VERB
    testflags += --nocapture --test-threads=1
endif

test:
	cargo test -- $(testflags) $(ARGS)

clippy:
	cargo clippy --workspace --tests
