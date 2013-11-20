#!/usr/bin/env ruby
require "#{File.dirname(__FILE__)}/fishers_exact_test.rb"
require "#{File.dirname(__FILE__)}/search_variants_util.rb"
FET = Rubystats::FishersExactTest.new

tumor_pileup = ARGV[0]
normal_pileup = ARGV[1]
out = ARGV[2]
@tumor_freq_th = 0.05
@normal_freq_th = 0.05
@tumor_depth_th = 8
@normal_depth_th = 8
@min_variant_num = 4
@base_qual_th = 15 + 33
@germline_freq_th = [0.25, 0.75]

go = open(out + "cand_hetero_germline", 'w')
so = open(out + "cand_somatic", 'w')

class NoMoreLinesInNormal < Exception
end

def get_new_record(f)
  l = f.gets
  raise NoMoreLinesInNormal if l.nil?
  Record.new(l)
end

open(tumor_pileup) { |tf|
  open(normal_pileup) { |nf|
    nr = nil
    count = 0
    prev_chr = nil
    begin
      nr = get_new_record(nf)
      while tl = tf.gets
        tr = Record.new(tl)
        count += 1
        if count % 100000 == 0
          warn [count, tr.chr, tr.pos].join("\t")
          GC.start
        end
        if tr.chr != nr.chr
          if prev_chr.nil?
            raise Exception, "We met different chrs in the first place."
          else
            if prev_chr != nr.chr # normal chr changed
              warn "err1: no matched locus in normal for #{tr.chr}:#{tr.pos}. skip!" if $DEBUG;next
            else
              while (nr.chr != tr.chr); nr = get_new_record(nf); end
            end
          end
        end
        while tr.pos > nr.pos # to catch the locus in tumor
          nr = get_new_record(nf)
          break if tr.chr != nr.chr # normal passed tumor
        end
        (warn "err2: no matched locus in normal for #{tr.chr}:#{tr.pos}. skip!" if $DEBUG;next) if tr.chr != nr.chr or tr.pos != nr.pos
        #print [tr.chr, tr.pos].join(",") + "|"
        prev_chr = tr.chr
        next if tr.ref == "N"
        next unless tr.depth >= @tumor_depth_th
        next unless nr.depth >= @normal_depth_th
        next if tr.bases.gsub(/(\^.)|[\.,\$]/,'').length < @min_variant_num
        warn "*" * 20 if $DEBUG
        warn tr if $DEBUG
        warn nr if $DEBUG
        tb = filter_and_bases(tr)
        nb = filter_and_bases(nr)
        tumor_depth = tb.inject(0) {|memo, o| o[0][0] == "-" ? memo : memo + o[1].sum }
        normal_depth = nb.inject(0) {|memo, o| o[0][0] == "-" ? memo : memo + o[1].sum }
        ref_tumor = tb.delete tr.ref
        ref_normal = nb.delete tr.ref
        if ref_tumor.nil? or ref_normal.nil?
          warn "warning: ref_tumor or ref_normal is nil. skip!"
          warn [tr.chr, tr.pos]
          warn [:tb, tb]
          warn [:nb, nb]
          next
        end
        tumor_depth = filtered_depth(tr)
        normal_depth = filtered_depth(nr)
        germ_out = []
        som_out = []
        for to in tb
          warn to if $DEBUG
          next if to[1].inject(:+) == 0
          oo = nb[to[0]] || [0, 0]
          if to[0][0] == "+"
            ref = "-"
            obs = to[0][1..-1]
            pos = tr.pos
            fin = tr.pos
          elsif to[0][0] == "-"
            ref = to[0][1..-1]
            obs = "-"
            pos = tr.pos
            fin = tr.pos + ref.length - 1
          else #base
            pos = tr.pos - 1
            fin = pos
            ref = tr.ref
            obs = to[0]
          end
          t_ref_c = ref_tumor.sum
          t_obs_c = to[1].sum
          n_ref_c = ref_normal.sum
          n_obs_c = oo.sum
          t_freq = t_obs_c.to_f / tumor_depth
          t_freq2 = t_obs_c.to_f / (t_obs_c.to_f + t_ref_c.to_f) # exclude other variants
          n_freq = (n_obs_c == 0) ? 0 : n_obs_c.to_f / normal_depth
          n_freq2 = (n_obs_c == 0) ? 0 : n_obs_c.to_f / (n_obs_c.to_f + n_ref_c.to_f) # exclude other variants
          t_st = to[1][0] / to[1].sum.to_f
          n_st = oo.sum == 0 ? "-" : (oo[0] / oo.sum.to_f)
          # filter
          warn [tr.chr, pos, fin, ref, obs, t_ref_c, t_obs_c, n_ref_c, n_obs_c, t_freq, t_st, n_freq, n_st].join("\t") if $DEBUG
          warn "check @min_variant_num" if $DEBUG
          next unless t_obs_c >= @min_variant_num
          if t_freq.between?(@germline_freq_th[0], @germline_freq_th[1]) and n_freq.between?(@germline_freq_th[0], @germline_freq_th[1])
            next unless t_freq2.between?(@germline_freq_th[0], @germline_freq_th[1]) and n_freq2.between?(@germline_freq_th[0], @germline_freq_th[1])
            next unless t_st.between?(0.1, 0.9) and n_st.between?(0.1, 0.9)
            germ_out.push [tr.chr, pos, fin, ref, obs, t_ref_c, t_obs_c, n_ref_c, n_obs_c, t_freq, t_st, n_freq, n_st].join("\t")
          elsif t_freq >= @tumor_freq_th and n_freq <= @normal_freq_th
            fs = FET.calculate(t_ref_c, t_obs_c, n_ref_c, n_obs_c)[:twotail]
            t_avg_bq = average_base_qualities(tr)
            n_avg_bq = average_base_qualities(nr)
            t_ref_bq = t_avg_bq[tr.ref]
            n_ref_bq = n_avg_bq[tr.ref]
            if ref != "-" and obs != "-"
              t_obs_bq = t_avg_bq[tr.ref]
              n_obs_bq = n_avg_bq[tr.ref]
            else
              t_obs_bq = "-"
              n_obs_bq = "-"
            end
            som_out.push [tr.chr, pos, fin, ref, obs, t_ref_c, t_obs_c, n_ref_c, n_obs_c, t_freq, t_st, n_freq, n_st, t_ref_bq, t_obs_bq, n_ref_bq, n_obs_bq, -Math.log10(fs)].join("\t")
          end
        end
        so.puts som_out.join("\n") if som_out.length == 1
        warn "multiple somatic variants in a locus: #{som_out};skip!" if som_out.length > 1
        go.puts germ_out[0] if germ_out.length == 1
        warn "multiple germline variants in a locus: #{germ_out};skip!" if germ_out.length > 1
      end
    rescue NoMoreLinesInNormal
      warn "hit the bottom of normal file. finish."
    end
  }
}
