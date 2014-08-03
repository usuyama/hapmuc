def rbeta(n, a, b)
  `R --vanilla --args #{n} #{a} #{b} rbeta < #{File.dirname(__FILE__)}/rbeta.R`
  open("rbeta").readline.chomp.split("\t").map {|x| x.to_f }
end

def gaussian(mu, sig) # N(mu, sig^2)
  t = -6.0
  12.times { t += rand() }
  t /= 4.0
  (t * sig) + mu
end

class Variant
  attr_reader :type, :ref_pos, :length, :str, :seq
  attr_accessor :hap_pos
  def initialize(pos, str)
    @str = str
    @ref_pos = pos
    @hap_pos = pos
    if /(\w)=>(\w)/ =~ str
      @type = :snp
      @seq = $2
      @length = 1
    elsif /\+(\w+)/ =~ str
      @type = :ins
      @seq = $1
      @length = $1.length
    elsif /-(\w+)/ =~ str
      @type = :del
      @seq = $1
      @length = $1.length
    end
  end
  def inspect
    "#{@ref_pos}, #{@hap_pos} " + @str
  end
end

class Haplotype
  attr_reader :variants, :seq, :punc
  def initialize(vars=[])
    @variants = vars.sort {|a,b| a.ref_pos <=> b.ref_pos }
    gen_seq
    gen_punc
  end

  def inspect
    @seq
  end

  def gen_punc
    @punc = []
    for v in @variants
      if v.type == :snp
        next
      elsif v.type == :ins
        @punc.push([:is, v.hap_pos])
        @punc.push([:ie, v.hap_pos+v.length])
      elsif v.type == :del
        @punc.push([:ds, v.hap_pos])
        @punc.push([:de, v.hap_pos+v.length])
      end
    end
  end

  def change_hap_pos(start, length)
    start.upto(@variants.size-1) do |i|
      @variants[i].hap_pos += length
    end
  end

  def seek_next(pos, rest)
    left = nil
    right = nil
    old = nil
    for x in @punc
      if x.last > pos
        right = x
        if !old.nil?
          left = old
        end
        break
      end
      old = x
    end
    if left.nil?
      if !right.nil?
        l = right.last - pos
        l = (l < rest) ? l : rest
        return [l, "#{l}M"]
      else
        [rest, "#{rest}M"]
      end
    end
    if right.nil?
      [rest, "#{rest}M"]
    else
      l = right.last - pos
      l = l > rest ? rest : l
      if  left.first == :ie || left.first == :de
        [l, "#{l}M"]
      elsif left.first == :is
        [l, "#{l}I"]
      else
        l = right.last - pos
        [l, "#{l}D"]
      end
    end
  end

  def gen_cigar(start)
    cp = start
    rest = READ_SIZE
    cigar = ""
    while(rest > 0)
      x = seek_next(cp, rest)
      cigar += x.last
      if x.last[-1,1]!="D"
        rest -= x.first
      end
      cp += x.first
    end
    cigar
  end


  def gen_seq
    @seq = REF.dup
    i = 0
    for v in @variants
      i += 1
      if v.type == :snp
        @seq[v.hap_pos-1] = v.seq
      elsif v.type == :ins
        @seq = @seq[0..v.hap_pos-1] + v.seq + @seq[v.hap_pos..@seq.size-1]
        change_hap_pos(i, v.length)
      elsif v.type == :del
        @seq = @seq[0..v.hap_pos-1] + @seq[v.hap_pos+v.length..@seq.size-1]
        change_hap_pos(i, -v.length)
      end
    end
  end

  def sanger_phred(pos)
    return 20
    high = 30
    offset = 20
    high_region = 40
    slope = (high-offset) / (100-40).to_f
    if pos < high_region
      high
    else
      offset + ((100-pos) * slope).to_i
    end
  end

  def phred_to_prob(phred_score)
    10 ** (-phred_score/10.0)
  end

  def phred_to_ascii(phred_score)
    (phred_score+33).chr
  end

  def gen_read(pos=nil)
    pos = rand(REF.length-READ_SIZE-1) if pos.nil?
    read = @seq[pos, READ_SIZE]
    qualities = ""
    0.upto(read.size-1) do |i|
      phred = sanger_phred(i)
      qualities += (phred_to_ascii(phred))
      if(rand() < phred_to_prob(phred))
        alt = read[i,1]
        while(alt==read[i,1])
          #alt = DNA[rand(4)]
          current = DNA.index(alt)
          alt = DNA[(current+1)%4] #決まった方向にエラーを起こす
        end
        read[i] = alt
      end
    end
    [read, qualities, gen_cigar(pos), pos]
  end

  def gen_paired_reads(pos)
    start = pos - READ_SIZE + rand(READ_SIZE)
    to_forward = rand() > 0.5 # direction of strand
    ins_size = gaussian(300, 100)
    if to_forward
      start1 = start
      start2 = start + ins_size - READ_SIZE
    else
      start2 = start
      start1 = start - ins_size + READ_SIZE
    end
    r1 = gen_read(start1)
    r2 = gen_read(start2)
    rev_r2 = [r2.first.gsub(/A/, "t").gsub(/T/, "a").gsub(/G/,"c").gsub(/C/,"g").upcase.reverse, r2[1].reverse, r2[2].scan(/\d+\w/).reverse.join, r2[3], r2[0]]
    [r1, rev_r2]
  end

  File.open(File.expand_path(File.dirname($0)) + "/header", "r") do |f|
    HEADER = f.read()
  end

  def fastq(seq_name, seq, qual)
    "@#{seq_name}\n#{seq}\n+#{seq_name}\n#{qual}\n"
  end

  def self.pair_to_sam(id, r)
    pos = r[0][3].to_i + POS + 2
    pair_pos = r[1][3].to_i + POS + 2
    ins = pair_pos - pos + 100
    out = "#{id}\t99\t#{CHR}\t#{pos}\t60\t100M\t=\t#{pair_pos}\t#{ins}\t#{r[0][0]}\t#{r[0][1]}\n"
    out += "#{id}\t147\t#{CHR}\t#{pair_pos}\t60\t100M\t=\t#{pos}\t#{-ins}\t#{r[1][4]}\t#{r[1][1]}\n"
    out
  end

  $count = 0
  def gen_paired_fastq(pos)
    $count += 1
    pair = gen_paired_reads(pos, pos)
    #[0,1].map{|i| fastq("#{$count}_#{pair[i][2]}/#{i+1}", pair[i][0], pair[i][1])}
    seq_name = "#{$count}_#{pair[0][2]}"
    [0,1].map{|i| fastq("#{seq_name}/#{i+1}", pair[i][0], pair[i][1])}
  end

end
