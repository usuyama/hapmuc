#!/usr/bin/env ruby
prev_symbol=nil
buf=[]

@min_distance = 0 # 0 then no filter

def to_symbol(r)
  if r[3] == "-"
    "#{r[1]}:+#{r[4]}"
  elsif r[4] == "-"
    "#{r[1]}:-#{r[3]}"
  else
    "#{r[1]}:#{r[3]}=>#{r[4]}"
  end
end

@cand_somatic_column_size = 20

def make_output(buf)
  if buf.length == 1 and buf[0].length < 31
    (buf[0].push "-").join("\t")
  else
    pos = buf[0][1].to_i
    out = buf[0][0..(@cand_somatic_column_size-1)]
    warn (pos.to_s + " " + ("*" * 10)) if $DEBUG
    warn out.join("\t") if $DEBUG
    out.push(buf.select {|x|
      g = x[@cand_somatic_column_size..-1]
      if (g[1].to_i - pos).abs < @min_distance
        warn "skip: too close germline var.: " + x.join("\t")
        false
      else;true;end
    }.map {|x|
      g = x[@cand_somatic_column_size..-1]
      to_symbol(g)
    }.join(','))
    out.join("\t")
  end
end

while l = gets
  r = l.chomp.split("\t")
  if prev_symbol != to_symbol(r)
    unless prev_symbol.nil?
      out = make_output(buf)
      puts out
    end
    buf = [r]
  else
    buf.push r
  end
  prev_symbol = to_symbol(r)
end

out = make_output(buf)
puts out
