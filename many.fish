#!/bin/fish
cargo build -r --examples

time for file in (find $argv[1] -name '*.xtc' | sort)
	echo '====>' $file && target/release/examples/compare $file || break
end
