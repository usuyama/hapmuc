require "#{File.dirname(__FILE__)}/seq_lib.rb"
refs = open("#{File.dirname(__FILE__)}/../tests/data/random_ref.fasta").read.split(">")[1..-1]#.split("\n")[1..-1].join("")
chr = rand(refs.size)
ref = refs[chr].split("\n")
REF = ref[1..-1].join("")
CHR = ref[0]
DNA = ["A", "T", "G", "C"]
READ_SIZE=100
POS = -1

depth = ARGV[0].to_i # roughly
tp = ARGV[1].to_f # tumor purity
@out_f = ARGV[2]
negative = !ARGV[3].nil? # generate a negative sample
no_germline_var = !ARGV[4].nil?
BASE_QUAL = 20

@fn = [0.5, 0.5, 0.0, 0.0] # h1, h2, h3, h4
@ft = [0.5 - tp/2, 0.5, tp/2, 0.0]

if negative
  @ft = [0.5 - tp/4, 0.5 - tp/4, tp/4, tp/4]
end
if no_germline_var
  @ft = [1.0 - tp/2, 0.0, tp/2, 0.0]
end

p [:normal_freq, @fn]
p [:tumor_freq, @ft]

def flip(pos)
  i = DNA.index(REF[pos])
  return i == DNA.length-1 ? DNA[0] : DNA[i+1]
end

snp_pos1 = 2000 + rand(1000)
snp_pos2 = 2000 + rand(1000)
snp_pos3 = 2000 + rand(1000)
snp1=Variant.new(snp_pos1+1, "#{REF[snp_pos1]}=>#{flip(snp_pos1)}")
snp2=Variant.new(snp_pos2+1, "#{REF[snp_pos2]}=>#{flip(snp_pos2)}")
snp3=Variant.new(snp_pos3+1, "#{REF[snp_pos3]}=>#{flip(snp_pos3)}")
@h1 = Haplotype.new([snp1])
@h2 = Haplotype.new([snp2, snp3])
sm_pos = snp_pos1 + rand(190) - 95
sm_pos += ((sm_pos - snp_pos1).abs < 5) ? 20 : 0
sm=Variant.new(sm_pos+1,"#{REF[sm_pos]}=>#{flip(sm_pos)}")
warn [snp1, snp2, snp3, sm].inspect
@h3 = Haplotype.new([sm, snp1])
@h4 = Haplotype.new([sm, snp2, snp3])
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

@normal_reads = []
index = 0
depth.to_i.times do
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
  pair = @haps[i].gen_paired_reads(sm_pos)
  normal.puts Haplotype.pair_to_sam(index, pair)
end

index = 0
@tumor_reads = []
depth.to_i.times do
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
  pair = @haps[i].gen_paired_reads(sm_pos)
  tumor.puts Haplotype.pair_to_sam(index, pair)
end

