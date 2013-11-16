class Record
  attr_accessor :chr, :pos, :ref, :depth, :bases, :quals
  def initialize(r)
    r = r.strip.split("\t")
    @chr = r[0]
    @pos = r[1].to_i
    @ref = r[2]
    @depth = r[3].to_i
    @bases = r[4]
    @quals = r[5]
  end
  def inspect;[chr, pos, ref, depth, bases, quals]; end
end

class CandVar
  attr_accessor :chr, :pos, :ref, :obs, :raw
  def initialize(r)
    @raw = r.chomp
    r = @raw.split("\t")
    @chr = r[0]
    @pos = r[1].to_i
    @ref = r[3]
    @obs = r[4]
  end
  def inspect;[chr, pos, ref, obs]; end
  def ==(v); v.inspect == self.inspect; end
  def format
    s = pos.to_s + ":"
    if ref == "-"
      s += "+#{obs}"
    elsif obs == "-"
      s += "-#{ref}"
    else
      s += "#{ref}=>#{obs}"
    end
    s
  end
end

class Array
  def sum
    self.inject(:+)
  end
end

def filtered_depth(r)
  r.quals.split("").map {|x| x.bytes.to_a[0] > @base_qual_th ? 1 : 0}.sum
end

def filter_and_bases(r)
  bases = r.bases.delete("$")
  quals = r.quals
  th = @base_qual_th
  i = 0
  data = {"A" => [0,0], "C" => [0,0], "G" => [0,0], "T" => [0,0]}
  while bases.length != 0
    warn [i, bases.length, bases] if $DEBUG
    f = bases[0]
    (bases = bases[2..-1];next) if f == "^"
    t = nil
    if bases =~ /^-(\d+)/
      t = "-" + bases[1+$1.length..$1.length+$1.to_i]
      bases = bases[$1.to_i+$1.length+1..-1]
      ut = t.upcase
      data[ut] = [0,0] if data[ut].nil?
      if t[-1] == ut[-1]
        data[ut][0] += 1
      else
        data[ut][1] += 1
      end
    elsif bases =~ /^\+(\d+)/
      t = "+" + bases[1+$1.length..$1.length+$1.to_i]
      bases = bases[$1.to_i+$1.length+1..-1]
      ut = t.upcase
      data[ut] = [0,0] if data[ut].nil?
      if t[-1] == ut[-1]
        data[ut][0] += 1
      else
        data[ut][1] += 1
      end
    else
      filter = (quals[i].bytes.to_a[0] > th)
      i += 1
      if f.upcase == "N" or f =="*" #n, N または del. なぜかdeleted baseにもqualがついてる
        bases = bases[1..-1]
      elsif f == "."
        bases = bases[1..-1]
        data[r.ref][0] += 1 if filter
      elsif f == ","
        bases = bases[1..-1]
        data[r.ref][1] += 1 if filter
      elsif f != "+" and f != "-"
        bases = bases[1..-1]
        if f == f.upcase
          data[f][0] += 1 if filter
        else
          data[f.upcase][1] += 1 if filter
        end
      end
    end
  end
  warn data if $DEBUG
  data
end

def average_base_qualities(r)
  bases = r.bases.delete("$")
  quals = r.quals
  th = @base_qual_th
  i = 0
  data = {"A" => [], "C" => [], "G" => [], "T" => []}
  while bases.length != 0
    warn [i, bases.length, bases] if $DEBUG
    f = bases[0]
    (bases = bases[2..-1];next) if f == "^"
    t = nil
    if bases =~ /^-(\d+)/
      t = "-" + bases[1+$1.length..$1.length+$1.to_i]
      bases = bases[$1.to_i+$1.length+1..-1]
      ut = t.upcase
      data[ut] = [0,0] if data[ut].nil?
      if t[-1] == ut[-1]
        data[ut][0] += 1
      else
        data[ut][1] += 1
      end
    elsif bases =~ /^\+(\d+)/
      t = "+" + bases[1+$1.length..$1.length+$1.to_i]
      bases = bases[$1.to_i+$1.length+1..-1]
      ut = t.upcase
      data[ut] = [0,0] if data[ut].nil?
      if t[-1] == ut[-1]
        data[ut][0] += 1
      else
        data[ut][1] += 1
      end
    else
      filter = (quals[i].bytes.to_a[0] > th)
      bq = quals[i].bytes.to_a[0] - 33
      i += 1
      if f.upcase == "N" or f =="*" #n, N または del. なぜかdeleted baseにもqualがついてる
        bases = bases[1..-1]
      elsif f == "."
        bases = bases[1..-1]
        #data[r.ref][0] += 1 if filter
        data[r.ref].push bq
      elsif f == ","
        bases = bases[1..-1]
        data[r.ref].push bq
        #data[r.ref][1] += 1 if filter
      elsif f != "+" and f != "-"
        bases = bases[1..-1]
        if f == f.upcase
          #data[f][0] += 1 if filter
          data[f].push bq
        else
          data[f.upcase].push bq
          #data[f.upcase][1] += 1 if filter
        end
      end
    end
  end
  warn data if $DEBUG
  out = {}
  "ATGC".split("").map {|x|
    bqs = data[x]
    if bqs.empty?
      out[x] = -1
    else
      out[x] = bqs.inject(:+).to_f / bqs.length
    end
  }
  return out
end
