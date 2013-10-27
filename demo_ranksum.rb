#!/usr/bin/env ruby

require_relative 'ranksum'

x = [8.50, 9.48, 8.65, 8.16, 8.83, 7.76, 8.63]
y = [8.27, 8.20, 8.25, 8.14, 9.00, 8.10, 7.20, 8.32, 7.70]
p = Ranksum.wilcoxon(x, y, nil)
puts "p = #{p}"

