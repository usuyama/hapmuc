#!/usr/bin/env ruby
# genenerate a random reference sequence
# 51 bases in 1 line
puts ">chr1"
bases = ['A', 'T', 'G', 'C']
100.times do
  51.times do
    print bases[rand(4)]
  end
  puts ''
end
