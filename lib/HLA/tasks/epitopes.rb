require 'rbbt/sources/organism'
Workflow.require_workflow "Sequence"
module HLA

  #input :fuseq_fusions, :tsv, "Fusions from FuSeq"
  #input :organism, :string, "Organism code", nil
  dep HTS, :FuSeq_process
  task :FuSeq_transcripts => :tsv do 
    fuseq_fusions = step(:FuSeq_process).load
    organism = recursive_inputs[:organism]

    transcript_sequence = Organism.transcript_sequence(organism).tsv :persist => true
    transcript_5utr = Organism.transcript_5utr(organism).tsv :persist => true
    transcript_3utr = Organism.transcript_3utr(organism).tsv :persist => true

    dumper = TSV::Dumper.new :key_field => "Fusion", :fields => ["CDNA", "Type", "Breakpoint", "5UTR", "3UTR"], :type => :double
    dumper.init
    TSV.traverse fuseq_fusions, :into => dumper do |fusion,values,fields|
      new = []
      Misc.zip_fields(values).each do |entry|
        info = Misc.zip2hash(fields, entry)

        strand5, strand3 = info.values_at "strand5", "strand3"

        next if strand5 == "-" and strand3 == "+"

        chr5, pos5 = info.values_at "chrom5", "brpos5"
        genomic_position5 = [chr5, pos5] * ":"
        transcript_offsets5 = Sequence.job(:transcript_offsets, nil, :positions => [genomic_position5], :organism => organism).run[genomic_position5]

        chr3, pos3 = info.values_at "chrom3", "brpos3"
        genomic_position3 = [chr3, pos3] * ":"
        transcript_offsets3 = Sequence.job(:transcript_offsets, nil, :positions => [genomic_position3], :organism => organism).run[genomic_position3]

        type = [strand5, strand3] * ":"
        case type
        when "+:+"
          new_start = transcript_offsets5.collect do |toffset|
            transcript, offset = toffset.split(":")
            offset = offset.to_i
            sequence = transcript_sequence[transcript]
            utr5 = transcript_5utr[transcript]
            next if utr5.nil?
            [sequence[0..offset-1], offset, utr5]
          end
          new_end = transcript_offsets3.collect do |toffset|
            transcript, offset = toffset.split(":")
            offset = offset.to_i
            sequence = transcript_sequence[transcript]
            utr3 = transcript_3utr[transcript]
            next if utr3.nil?
            [sequence[offset..-1], utr3]
          end
          new_start.compact.each do |ns,os,utr5|
            new_end.compact.each do |ne,utr3|
              seq = [ns, ne] * ""
              new << [seq, type, os, utr5, utr3]
            end
          end
        when "-:-"
          new_start = transcript_offsets5.collect do |toffset|
            transcript, offset = toffset.split(":")
            offset = offset.to_i
            sequence = transcript_sequence[transcript]
            utr5 = transcript_5utr[transcript]
            next if utr5.nil?
            [sequence[0..offset-1], offset, utr5]
          end
          new_end = transcript_offsets3.collect do |toffset|
            transcript, offset = toffset.split(":")
            offset = offset.to_i
            sequence = transcript_sequence[transcript]
            utr3 = transcript_3utr[transcript]
            next if utr3.nil?
            [sequence[offset..-1], utr3]
          end
          new_start.compact.each do |ns,os,utr5|
            new_end.compact.each do |ne,utr3|
              seq = [ns, ne] * ""
              new << [seq, type, os, utr5, utr3]
            end
          end
        when "+:-"
        end
      end
      next if new.empty?
      new.extend MultipleResult
      [fusion, Misc.zip_fields(new)]
    end
  end

  dep :FuSeq_transcripts
  task :FuSeq_transcript_evidence => :tsv do
  end

  input :mutated_isoforms, :array, "Mutated Isoforms"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :mutated_isoform_sequences => :tsv do |mis,organism|
    mi_sequence = Organism.protein_sequence(organism).tsv :persist => true
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Wildtype", "Pos", "WT", "Alt"], :type => :list, :organism => organism
    dumper.init
    TSV.traverse mis, :type => :array, :into => dumper do |mi|
      protein, _sep, change = mi.partition(":")
      wt, pos, alt = change.partition /\d+/
      wt_sequence = mi_sequence[protein]
      raise "No sequence for protein: #{ protein }" if wt_sequence.nil?
      raise "Wildtype aa at position #{pos} is #{wt_sequence[pos.to_i - 1]} not #{wt}" if wt != wt_sequence[pos.to_i - 1] 
      [mi, [wt_sequence, pos, wt, alt]]
    end
  end

  dep :mutated_isoform_sequences
  input :flank_size, :integer, "Flank size", 20
  task :mutation_flanking_sequence => :tsv do |flank_size|
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Lost", "Position of Change", "N Flank", "Change", "C Flank"], :type => :list
    dumper.init
    TSV.traverse step(:mutated_isoform_sequences), :into => dumper do |mi,values|
      mi = mi.first if Array === mi

      wt_sequence, pos, wt, alt = values
      pos = pos.to_i

      mut_sequence = wt_sequence.dup
      next if alt == "*"

      aa_size = 1
      lost = ""
      if m = alt.match(/(Indel|FrameShift)\((.*)\)/i)
        type = m[1]
        aas = m[2]

        next if aas.nil?
        next if aas == ""
        next if aas == "?"

        aas = "-" + aas unless %(+ -).include? aas[0]

        aas_a = aas.split("")
        while aas_a.first == '-'
          lost << mut_sequence[pos - 1]
          mut_sequence[pos - 1] = ""
          aas_a.shift
        end

        if aas_a.first == '+'
          pos += 1
          aas_a.shift
        end

        if type == "Indel"
          aa_size = aas_a.length
          mut_sequence[pos - 1] = (aas_a * "") + mut_sequence[pos - 1]  if aas_a
        else
          aa_size = 0
          lost = ""
          mut_sequence[pos - 1..-1] = aas_a * ""
        end
      else
        type = 'SNV'
        lost << mut_sequence[pos - 1]
        mut_sequence[pos - 1] = alt
      end

      n_start = pos-flank_size-1
      n_start = 0 if n_start < 0

      if pos > 1
        n_flank = mut_sequence[n_start..pos-2]
      else
        n_flank = ""
      end

      if type === "Indel" || aa_size > 0 # Indel
        pos_end = pos + aa_size - 1
        pos_str = mut_sequence[pos-1..pos_end-1]
        c_flank = mut_sequence[pos_end..pos+flank_size+aa_size-2]
      else # FrameShift
        pos_str = mut_sequence[pos-1..-1]
        pos_str << "-" unless pos_str[-1] == "*"
        c_flank = wt_sequence[(pos-1)..pos+flank_size-2]
      end

      iii mi
      iif [flank_size, n_flank, + n_flank]
      n_flank = "-" * (flank_size - n_flank.length) + n_flank
      c_flank = c_flank + "-" * (flank_size - c_flank.length)
      [mi, [lost, pos, n_flank, pos_str, c_flank]]
    end
  end

  dep :mutation_flanking_sequence, :compute => :produce
  input :sizes, :array, "Peptide size", [9]
  task :epitopes => :tsv do |sizes|
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Wildtype", "Mutated",  "Offset"], :type => :double
    dumper.init
    TSV.traverse step(:mutation_flanking_sequence), :into => dumper do |mi,info|
      mi = mi.first if Array === mi
      lost, change_position, nterm, change, cterm = info

      frameshift = mi.include? "FrameShift"

      offset = nterm.length
      mut_peptide = nterm + change + cterm
      wt_peptide = nterm + lost + cterm

      epitopes = []
      offset = nterm.length
      diff = change.length - lost.length

      sizes.each do |size|
        size = size.to_i
        middle = change.length + (size - change.length) / 2 
        steps = size + change.length - 1
        steps.times do |step|
          index = offset + step - (size - 1)
          mut_epitope = mut_peptide[index..index+size - 1]
          next if mut_epitope.include? '-'
          next if mut_epitope.include? '*'

          left_wt_epitope = wt_peptide[index..index+size - 1]
          left_wt_epitope = "" if left_wt_epitope.nil? || left_wt_epitope.nil?
          right_wt_epitope = wt_peptide[index-diff..index+size-diff-1]
          right_wt_epitope = "" if frameshift || right_wt_epitope.nil?

          right_wt_epitope = "" if right_wt_epitope.include? "*"
          left_wt_epitope = "" if left_wt_epitope.include? "*"
          right_wt_epitope = "" if right_wt_epitope.length != size
          left_wt_epitope = "" if left_wt_epitope.length != size

          left_wt_epitope_common = 0
          size.times do |i|
            left_wt_epitope_common += 1 if left_wt_epitope[i] == mut_epitope[i]
          end

          right_wt_epitope_common = 0
          size.times do |i|
            right_wt_epitope_common += 1 if right_wt_epitope[i] == mut_epitope[i]
          end

          wt_epitope = if right_wt_epitope_common > left_wt_epitope_common
                         right_wt_epitope
                       else
                         left_wt_epitope
                       end

          wt_epitope = "" if [right_wt_epitope_common, left_wt_epitope_common].max <= [3, size / 2].max

          next if mut_epitope == wt_epitope

          epitopes << [wt_epitope, mut_epitope, step]
        end
      end

      [mi, Misc.zip_fields(epitopes)]
    end
  end

  dep :mutated_isoform_sequences
  input :sizes, :array, "Peptide size", [9]
  task :epitopes_old => :tsv do |sizes|
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Wildtype", "Mutated", "Offset", "Flank"], :type => :double
    dumper.init
    TSV.traverse step(:mutated_isoform_sequences), :into => dumper do |mi,values|
      mi = mi.first if Array === mi
      epitopes = []
      wt_sequence, pos, wt, alt = values
      mut_sequence = wt_sequence.dup
      next if alt == "*"

      aa_size = 1
      if m = alt.match(/(Indel|FrameShift)\((.*)\)/i)
        type = m[1]
        aas = m[2]
        aas_a = aas.split("")
        while aas_a.first == '-'
          mut_sequence[pos.to_i - 1] = ""
          aas_a.shift
        end
        if type == "Indel"
          aa_size = aas_a.length
          mut_sequence[pos.to_i - 1] = (aas_a * "") + mut_sequence[pos.to_i - 1]  if aas_a
        else
          mut_sequence[pos.to_i - 1..-1] = aas_a * ""
        end
      else
        type = 'SNV'
        mut_sequence[pos.to_i - 1] = alt
      end

      sizes.each do |size|
        size = size.to_i
        extra = (aa_size - 1)
        (size + extra).times do |i|
          start = pos.to_i - 1 - i + (aa_size - 1)
          next if start < 0

          eend = start + size - 1
          mut_epitope = mut_sequence[start..eend]
          mut_epitope = mut_epitope[0..-2] if mut_epitope[-1] == "*"
          next if mut_epitope.length < size
          next if mut_epitope.include?"*"

          wt_epitope1 = wt_sequence[start..eend]
          wt_epitope2 = wt_sequence[(start-extra)..(eend-extra)]
          wt_epitope = (i + 1) > ((size + extra)/2).round ?  wt_epitope1 : wt_epitope2
          wt_epitope = wt_epitope2 if wt_epitope.nil?

          wt_epitope = wt_epitope[0..-2] if wt_epitope && wt_epitope[-1] == "*"
          wt_epitope = nil if wt_epitope && wt_epitope.length < size
          wt_epitope = nil if wt_epitope && wt_epitope.include?("*")

          next if wt_epitope == mut_epitope
          epitopes << [wt_epitope, mut_epitope, i]
        end

        if type == "FrameShift"
          remain = mut_sequence.length - pos.to_i - size + 1
          remain.times do |i|
            start = pos.to_i + i
            eend = start + size - 1
            mut_epitope = mut_sequence[start..eend]
            mut_epitope = mut_epitope[0..-2] if mut_epitope[-1] == "*"
            next if mut_epitope.length < size
            next if mut_epitope.include?("*")
            epitopes << [nil, mut_epitope, start + 1]
          end
        end
      end


      next if epitopes.empty?
      [mi, Misc.zip_fields(epitopes)]
    end
  end

  dep :epitopes
  dep :epitopes_old
  task :cmp_epi => :array do
    tsv1 = step(:epitopes).load
    tsv2 = step(:epitopes_old).load
    diff = []
    tsv1.each do |k,v1|
      v2 = tsv2[k]
      p1 = v1[0].zip(v1[1]).sort_by{|a,b| b}
      p2 = v2[0].zip(v2[1]).sort_by{|a,b| b}
      if p1 != p2
        eee "diff" 
        wwww k
        wwww p1
        wwww p2
        wwww p1 - p2
        wwww p2 - p1
        diff << k
      end
    end
    diff
  end

  dep :epitopes
  input :alleles, :array, "List of alleles"
  task :mhcFlurry => :tsv do |alleles|
    str = 'allele,peptide' << "\n"
    TSV.traverse step(:epitopes), :into => str do |mi,values|
      mi = mi.first if Array === mi
      wt_epitopes, mut_epitopes, offsets = values
      next if mut_epitopes.nil? or mut_epitopes.compact.empty?
      all_epitopes = (wt_epitopes + mut_epitopes).compact.uniq
      all_epitopes.reject!{|e| e.empty?}

      res = []

      res = alleles.collect{|a| all_epitopes.compact.uniq.collect{|epitope| [a,epitope] * "," } }.flatten

      res * "\n" + "\n"
    end

    input = file('input.csv')
    output = file('output.csv')

    cpus = config("cpus", :mchflurry, :neo_epitopes, :default => 3)
    Open.write(input, str)
    cmd = "env MHCFLURRY_DOWNLOADS_DIR='#{Rbbt.share.databases.mhcflurry.find}' mhcflurry-predict"
    cmd << " --threads #{cpus}" if cpus
    CMD.cmd_log(cmd + " #{ input } --out #{output}" )

    mhcflurry_scores = {}
    TSV.traverse output, :type => :array do |line|
      next if line =~ /mhcflurry_prediction/
      allele, epitope, pred, low, high, percent = line.split(",")
      mhcflurry_scores[[allele, epitope]] = [pred, low, high, percent]
    end

    dumper = TSV::Dumper.new :key_field => "Mutation ID", :type => :list, :fields => ["Mutated Isoform", "Allele", "Mutated epitope", "Wildtype epitope", "Offset",
                                                                      "Mutated mhcflurry prediction", "Mutated mhcflurry low", "Mutated mhcflurry high", "Mutated mhcflurry percentile", 
                                                                      "Wildtype mhcflurry prediction", "Wildtype mhcflurry low", "Wildtype mhcflurry high", "Wildtype mhcflurry percentile"]
    dumper.init
    TSV.traverse step(:epitopes), :into => dumper do |mi,values|
      mi = mi.first if Array === mi
      wt_epitopes, mut_epitopes, offsets = values
      next if mut_epitopes.nil? or mut_epitopes.compact.empty?

      res = []
      Misc.zip_fields(values).each do |wt_epitope,mut_epitope,offset|
        alleles.each do |allele|
          size = mut_epitope.length
          id = [mi, size, offset, allele] * "."
          wt_scores = mhcflurry_scores[[allele, wt_epitope]]
          mut_scores = mhcflurry_scores[[allele, mut_epitope]]
          wt_scores = [nil] * mut_scores.length if wt_scores.nil?
          res << [id, [mi, allele, mut_epitope, wt_epitope, offset] + mut_scores + wt_scores]
        end
      end
      res.extend MultipleResult
      res 
    end
  end

end
