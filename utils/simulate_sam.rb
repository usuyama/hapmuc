require "#{File.dirname(__FILE__)}/seq_lib.rb"
REF  = open("#{File.dirname(__FILE__)}/../tests/random_ref.fasta").read.split("\n")[1..-1].join("")
DNA = ["A", "T", "G", "C"]
READ_SIZE=100
CHR = "chr1"
POS = -1

depth = ARGV[0].to_i # roughly
tp = ARGV[1].to_f # tumor purity
@out_f = ARGV[2]
negative = !ARGV[3].nil? # generate a negative sample
BASE_QUAL = 20

@fn = [0.5, 0.5, 0.0, 0.0] # h1, h2, h3, h4
@ft = [0.5 - tp/2, 0.5, tp/2, 0.0]

if negative
  @ft = [0.5 - tp/4, 0.5 - tp/4, tp/4, tp/4]
end

p [:normal_freq, @fn]
p [:tumor_freq, @ft]

def flip(pos)
  i = DNA.index(REF[pos])
  return i == DNA.length-1 ? DNA[0] : DNA[i+1]
end

snp_pos = 1000 + rand(2000)
@h1 = Haplotype.new
snp=Variant.new(snp_pos+1, "#{REF[snp_pos]}=>#{flip(snp_pos)}")
@h2 = Haplotype.new([snp])
sm_pos = snp_pos + rand(190) - 95
sm_pos += ((sm_pos - snp_pos).abs < 5) ? 20 : 0
sm=Variant.new(sm_pos+1,"#{REF[sm_pos]}=>#{flip(sm_pos)}")
@h3 = Haplotype.new([sm])
@h4 = Haplotype.new([sm, snp])
`echo "#{POS + sm_pos + 1}" > ans` if tp > 0.0

if rand() > 0.5
  @haps = [@h1, @h2, @h3, @h4]
else
  @haps = [@h2, @h1, @h4, @h3]
end

tumor = File.open(@out_f + "tumor.sam", "w")
normal = File.open(@out_f + "normal.sam", "w")
tumor.puts Haplotype::HEADER
normal.puts Haplotype::HEADER

if sm_pos < snp_pos
  l_pos = sm_pos;r_pos = snp_pos;
else
  l_pos = snp_pos;r_pos = sm_pos;
end

@normal_reads = []
index = 0
(depth*(200 + r_pos-l_pos).to_f/READ_SIZE).to_i.times do
  r = rand()
  i = 0
  th = 0.0
  0.upto(@haps.size()-1) do
    th += @fn[i]
    if th >= r
      break
    else
      i += 1
    end
  end
  index+=1
  pair = @haps[i].gen_paired_reads(l_pos, r_pos)
  normal.puts Haplotype.pair_to_sam(index, pair)
end

index = 0
@tumor_reads = []
(depth*(200 + r_pos-l_pos).to_f/READ_SIZE).to_i.times do
  r = rand()
  i = 0
  th = 0.0
  0.upto(@haps.size()-1) do
    th += @ft[i]
    if th >= r
      break
    else
      i += 1
    end
  end
  index+=1
  pair = @haps[i].gen_paired_reads(l_pos, r_pos)
  tumor.puts Haplotype.pair_to_sam(index, pair)
end

