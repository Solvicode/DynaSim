.PHONY: test

test: 
	find . -name "test_*.zig" -exec zig test {} +

