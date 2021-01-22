require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/sources/organism'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/TFCheckpoint'

module TFCheckpoint
  extend Workflow


  input :organism, :string, "Organism code", "Hsa"
  input :remove_unknown, :boolean, "Remove genes that cannot be identified", false
  task :join => :tsv do |organism, remove_unknown|
    join = TSV.setup({}, :key_field => "Associated Gene Name", :fields => [], :type => :double)
    TFCheckpoint.root.glob("*").each do |file|
      tsv = TSV.open file

      fields = tsv.fields

      db = File.basename(file)
      tsv.fields = [db + ": ID used (#{fields.first})"] + fields[1..-1].collect{|f| [db, f] * ": "}
      next unless tsv.namespace.include? organism

      tsv = tsv.change_key "Associated Gene Name", :identifiers => Organism.identifiers(organism) unless tsv.key_field == "Associated Gene Name"

      join = join.attach tsv, :complete => true, :fields => tsv.fields, :one2one => false
    end

    name2entrez = Organism.identifiers(organism).index :target => "Entrez Gene ID", :persist => true
    join = join.add_field "Entrez Gene ID" do |k,values|
      name2entrez[k]
    end

    name2uni = Organism.identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
    join = join.add_field "UniProt/SwissProt Accession" do |k,values|
      name2uni[k]
    end

    name2ensembl = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :persist => true
    join = join.add_field "Ensembl Gene ID" do |k,values|
      name2ensembl[k]
    end
    join.namespace = organism
    
    join = join.select("Ensembl Gene ID"){|i| i.first } if remove_unknown

    join
  end

  dep :join, :organism => "Hsa", :jobname => "Hsa"
  dep :join, :organism => "Mmu", :jobname => "Mmu"
  dep :join, :organism => "Rno", :jobname => "Rno"
  task :all_orgs => :tsv do
    hsa, rno, mmu = dependencies.collect do |dep| 
      tsv = dep.load
      organism = tsv.namespace
      tsv.add_field "Organism" do
        organism
      end
      tsv.add_field "Entrez Taxa ID" do
        Organism.entrez_taxids(organism).produce.read
      end
    end
    
    hsa = hsa.attach rno, :complete  => true, :merge => true, :fields => rno.fields
    hsa = hsa.attach mmu, :complete  => true, :merge => true, :fields => mmu.fields
    hsa
  end

  dep :all_orgs
  task :ortho => :tsv do
    tsv = step(:all_orgs).load.to_double
    entrez_index = tsv.index :fields => ["Entrez Gene ID"]

    taxas = tsv.column("Entrez Taxa ID").values.flatten.uniq.compact.sort
    orthologs = {}
    TSV.traverse Rbbt.data["gene_orthologs"], :type => :array, :bar => true do |line|
      taxid, gene, rel, otaxid, ogene = line.split("\t")
      next unless taxas.include?(taxid) && taxas.include?(otaxid) && [taxid, otaxid].include?('9606')
      if taxid == '9606'
        taxid, otaxid = otaxid, taxid
        gene, ogene = ogene, gene
      end
      source = [taxid, gene] * ":"
      target = [otaxid, ogene] * ":"
      orthologs[source] = []
      orthologs[source] << target
    end
    
    new = tsv.annotate({})

    TSV.traverse tsv, :into => new do |name, values|
      organism, taxa, entrez = values.values_at "Organism", "Entrez Taxa ID", "Entrez Gene ID"
      ortho_key = [taxa, entrez] * ":"
      ovalues = orthologs[ortho_key]
      res = []
      res.extend MultipleResult
      if taxa != '9606' && ovalues 
        ovalues.each do |k|
          t, e = k.split(":")
          name = entrez_index[e]
          next if name.nil?
          res << [name, values]
        end
      else
        res << [name, values]
      end
      res
    end
    
  end

  dep :all_orgs
  task :ortho2 => :tsv do
    tsv = step(:all_orgs).load.to_double
    ensembl_index = tsv.index :fields => ["Ensembl Gene ID"]

    new = tsv.annotate({})

    organism_orthologs_all = {}
    organism_orthologs_all["Hsa"] = {}

    TSV.traverse tsv, :into => new do |name, values|
      organism, ensembl = values.values_at "Organism", "Ensembl Gene ID"
      organism = organism.first if Array === organism

      organism_orthologs = organism_orthologs_all[organism] ||= Organism.ortholog_Hsa(organism).tsv
      orthologs = organism_orthologs.values_at(*ensembl).compact.flatten
      res = []
      res.extend MultipleResult
      if ! organism.include?("Hsa") && orthologs.any?
        orthologs.each do |e|
          name = ensembl_index[e]
          next if name.nil?
          res << [name, values]
        end
      else
        res << [name, values]
      end

      res
    end
    
  end

  dep :all_orgs
  task :ortho3 => :tsv do
    tsv = step(:all_orgs).load.to_double

    ensembl_index = tsv.index :fields => ["Ensembl Gene ID"]
    entrez_index = tsv.index :fields => ["Entrez Gene ID"]

    new = tsv.annotate({})

    taxas = tsv.column("Entrez Taxa ID").values.flatten.uniq.compact.sort
    ortholog_entrez = {}
    TSV.traverse Rbbt.data["gene_orthologs"], :type => :array, :bar => true do |line|
      taxid, gene, rel, otaxid, ogene = line.split("\t")
      next unless taxas.include?(taxid) && taxas.include?(otaxid) && [taxid, otaxid].include?('9606')
      if taxid == '9606'
        taxid, otaxid = otaxid, taxid
        gene, ogene = ogene, gene
      end
      source = [taxid, gene] * ":"
      target = [otaxid, ogene] * ":"
      ortholog_entrez[source] = []
      ortholog_entrez[source] << target
    end
    organism_orthologs_all = {}
    organism_orthologs_all["Hsa"] = {}

    TSV.traverse tsv, :into => new do |name, values|
      organism, ensembl, taxa, entrez  = values.values_at "Organism", "Ensembl Gene ID", "Entrez Taxa ID", "Entrez Gene ID"
      organism = organism.first if Array === organism

      organism_orthologs = organism_orthologs_all[organism] ||= Organism.ortholog_Hsa(organism).tsv
      orthologs = organism_orthologs.values_at(*ensembl).compact.flatten

      ortho_keys = Misc.zip_fields([taxa, entrez]).collect{|p| p * ":"}

      ovalues = ortholog_entrez.values_at(*ortho_keys).flatten.uniq.compact


      res = []
      res.extend MultipleResult
      if ! organism.include?("Hsa") && orthologs.any?
        orthologs.each do |e|
          name = ensembl_index[e]
          next if name.nil?
          res << [name, values]
        end
      elsif taxa != '9606' && ovalues.any? 
        ovalues.each do |k|
          t, e = k.split(":")
          name = entrez_index[e]
          next if name.nil?
          res << [name, values]
        end
      else
        res << [name, values]
      end

      res
    end
    
  end
  
  dep :all_orgs
  task :ortho4 => :tsv do
    tsv = step(:all_orgs).load.to_double

    ensembl_index = tsv.index :fields => ["Ensembl Gene ID"]
    entrez_index = tsv.index :fields => ["Entrez Gene ID"]

    new = tsv.annotate({})

    taxas = tsv.column("Entrez Taxa ID").values.flatten.uniq.compact.sort
    ortholog_entrez = {}
    TSV.traverse Rbbt.data["gene_orthologs"], :type => :array, :bar => true do |line|
      taxid, gene, rel, otaxid, ogene = line.split("\t")
      next unless taxas.include?(taxid) && taxas.include?(otaxid) && [taxid, otaxid].include?('9606')
      if taxid == '9606'
        taxid, otaxid = otaxid, taxid
        gene, ogene = ogene, gene
      end
      source = [taxid, gene] * ":"
      target = [otaxid, ogene] * ":"
      ortholog_entrez[source] = []
      ortholog_entrez[source] << target
    end
    organism_orthologs_all = {}
    organism_orthologs_all["Hsa"] = {}

    TSV.traverse tsv, :into => new do |name, values|
      organism, ensembl, taxa, entrez  = values.values_at "Organism", "Ensembl Gene ID", "Entrez Taxa ID", "Entrez Gene ID"
      organism = organism.first if Array === organism

      organism_orthologs = organism_orthologs_all[organism] ||= Organism.ortholog_Hsa(organism).tsv
      orthologs = organism_orthologs.values_at(*ensembl).compact.flatten

      ortho_keys = Misc.zip_fields([taxa, entrez]).collect{|p| p * ":"}

      ovalues = ortholog_entrez.values_at(*ortho_keys).flatten.uniq.compact


      res = []
      res_ensembl = []
      res_entrez = []
      if ! organism.include?("Hsa") && orthologs.any?
        orthologs.each do |e|
          name = ensembl_index[e]
          next if name.nil?
          res_ensembl << [name, values]
        end
      elsif taxa != '9606' && ovalues.any? 
        ovalues.each do |k|
          t, e = k.split(":")
          name = entrez_index[e]
          next if name.nil?
          res_entrez << [name, values]
        end
      else
        res << [name, values]
      end

      if res.empty?
        if res_entrez.empty?
          res = res_ensembl
        elsif res_ensembl.empty?
          res = res_entrez
        else
          if res_ensembl.first.downcase == name.downcase
            res = res_ensembl
          else
            res = res_entrez
          end
        end
      end

      res.extend MultipleResult
      res
    end
    
  end


  dep :all_orgs
  task :ortho5 => :tsv do
    tsv = step(:all_orgs).load.to_double
    entrez_index = tsv.index :fields => ["Entrez Gene ID"]

    taxas = tsv.column("Entrez Taxa ID").values.flatten.uniq.compact.sort
    orthologs = {}

    mmu2hsa = Rbbt.data["human_mouse_orthologs.tsv"].tsv(:header_hash => "", :merge => true, :type => :flat, :key_field => 'mouse')
    rno2hsa = Rbbt.data["human_rat_orthologs.tsv"].tsv(:header_hash => "", :merge => true, :type => :flat, :key_field => 'rat')
    rno2mmu = Rbbt.data["mouse_rat_orthologs.tsv"].tsv(:header_hash => "", :merge => true, :type => :flat, :key_field => 'rat')

    i = 0
    new = tsv.annotate({})
    field_pos = ["Organism", "Entrez Taxa ID", "Entrez Gene ID"].collect{|f| tsv.fields.index(f) }
    TSV.traverse tsv, :into => new do |name, values|
      res = []
      res.extend MultipleResult

      used_names = []
      organismsl = values["Organism"]
      entrezl = values["Entrez Gene ID"]
      Misc.zip_fields([organismsl, entrezl]).each do |organism,entrez|
        next if organism.nil?

        case organism
        when 'Mmu'
          oe = (mmu2hsa[entrez] || []) - entrezl
          onames = entrez_index.values_at(*oe)
          if onames.any?
            onames.each do |oname|
              next if used_names.include? oname
              used_names << oname
              res << [oname, values]
            end
          else
            res << [name, values]
          end
        when 'Rno'
          oe = (rno2hsa[entrez] || []) - entrezl

          onames = entrez_index.values_at(*oe)

          if onames.any?
            onames.each do |oname|
              next if used_names.include? oname
              used_names << oname
              res << [oname, values]
            end
          else
            oe = (rno2mmu[entrez] || []) - entrezl
            onames = entrez_index.values_at(*oe)
            if onames.any?
              onames.each do |oname|
                next if used_names.include? oname
                used_names << oname
                res << [oname, values]
              end
            else
              res << [name, values]
            end
          end
        else
          res << [name, values]
        end
      end
      res.uniq!
      res
    end
    
  end

  dep :all_orgs
  task :ortho6 => :tsv do
    tsv = step(:all_orgs).load.to_double
    entrez_index = tsv.index :fields => ["Entrez Gene ID"]

    taxas = tsv.column("Entrez Taxa ID").values.flatten.uniq.compact.sort
    orthologs = {}

    mmu2hsa = Rbbt.data["human_mouse_orthologs.tsv"].tsv(:header_hash => "", :merge => true, :type => :flat, :key_field => 'mouse')
    rno2hsa = Rbbt.data["human_rat_orthologs.tsv"].tsv(:header_hash => "", :merge => true, :type => :flat, :key_field => 'rat')
    rno2mmu = Rbbt.data["mouse_rat_orthologs.tsv"].tsv(:header_hash => "", :merge => true, :type => :flat, :key_field => 'rat')
    
    manual_orthologs = TSV.xlsx(Rbbt.data["TFcheckpoint_orthologs_manual.xlsx"].find)

    manual_orthologs.process "mouse NCBI_D" do |v|
      v.first.to_i.to_s
    end

    manual_orthologs.process "rat NCBI_ID" do |v|
      v.first.to_i.to_s
    end

    manual_orthologs_mmu = manual_orthologs.index :target => 'Ortholog', :fields => ["mouse NCBI_D"]
    manual_orthologs_rno = manual_orthologs.index :target => 'Ortholog', :fields => ["rat NCBI_ID"]


    i = 0
    new = tsv.annotate({})
    field_pos = ["Organism", "Entrez Taxa ID", "Entrez Gene ID"].collect{|f| tsv.fields.index(f) }
    TSV.traverse tsv, :into => new do |name, values|
      res = []
      res.extend MultipleResult

      used_names = []
      organismsl = values["Organism"]
      entrezl = values["Entrez Gene ID"]
      Misc.zip_fields([organismsl, entrezl]).each do |organism,entrez|
        next if organism.nil?

        case organism
        when 'Mmu'
          oe = (mmu2hsa[entrez] || []) - entrezl
          onames = entrez_index.values_at(*oe)
          new_oname = manual_orthologs_mmu[entrez]
          onames << new_oname if new_oname
          if onames.any?
            onames.each do |oname|
              next if used_names.include? oname
              used_names << oname
              res << [oname, values]
            end
          else
            res << [name, values]
          end
        when 'Rno'
          oe = (rno2hsa[entrez] || []) - entrezl

          onames = entrez_index.values_at(*oe)
          new_oname = manual_orthologs_rno[entrez]
          onames << new_oname if new_oname

          if onames.any?
            onames.each do |oname|
              next if used_names.include? oname
              used_names << oname
              res << [oname, values]
            end
          else
            oe = (rno2mmu[entrez] || []) - entrezl
            onames = entrez_index.values_at(*oe)
            if onames.any?
              onames.each do |oname|
                next if used_names.include? oname
                used_names << oname
                res << [oname, values]
              end
            else
              res << [name, values]
            end
          end
        else
          res << [name, values]
        end
      end
      res.uniq!
      res
    end
    
  end

  dep :ortho6
  task :tidy => :tsv do
    tsv = dependencies.first.load
    fields = tsv.fields

    new_fields =<<-EOF.split("\n")
Entrez Taxa ID
Entrez Gene ID
UniProt/SwissProt Accession
Ensembl Gene ID
    EOF

    go_fields = {}
    fields.each do |f|
      go, num, org = f.split("_")
      next unless go == "go"
      go_fields[num] ||= []
      go_fields[num] << f
    end

    go_fields.each do |num,fields|
      term = fields.select{|f| f.include? "TERM"}
      evidence = fields.select{|f| f.include? "EVIDENCE"}

      tsv = tsv.add_field "GO #{num} terms" do |k,values|
        values.values_at(*term).flatten.compact
      end

      tsv = tsv.add_field "GO #{num} evidence" do |k,values|
        values.values_at(*evidence).flatten.compact
      end
    end

    fields = tsv.fields


    ids = fields.select{|f| f.include? "ID used"}
    go_fields = fields.select{|f| f=~/^GO /}
    old_go_fields = fields.select{|f| f=~/^go_/}

    features = fields.select{|f| f.include? ":"} - ids - go_fields - old_go_fields
    rest = fields - new_fields - ids - features - go_fields - old_go_fields

    tsv = tsv.reorder :key, new_fields + features + go_fields + ids 

    ids.sort.each do |field|
      db = field.split(": ").first
      next if db =~ /^go_/
      tsv.add_field "#{db} present" do |k,values|
        if values[field].any?
          [db]
        end
      end
    end

    tsv = tsv.reorder :key, tsv.fields - ids

    # Remove GO duplicates


    tsv.through do |k,values|
      %w(140223 3700 3712 43565 6355 6357 981).each do |goid|
        fterm = "GO #{goid} terms"
        fevidence = "GO #{goid} evidence"

        terms = values[fterm]
        evidence = values[fevidence]

        uniq = Misc.zip_fields([terms, evidence]).uniq
        terms = uniq.collect{|p| p[0] }
        evidence = uniq.collect{|p| p[1] }

        values[fterm] = terms
        values[fevidence] = evidence
      end

      id_fields = ["Entrez Taxa ID", "Entrez Gene ID", "UniProt/SwissProt Accession", "Ensembl Gene ID"]
      id_values = values.values_at *id_fields
      uniq_id_values = Misc.zip_fields(Misc.zip_fields(id_values).uniq)
      id_fields.zip(uniq_id_values).each do |field, value|
        values[field] = value
      end
    end

    tsv
  end

  dep :tidy
  task :filter_go => :tsv do
    tsv = step(:tidy).load
    present_fields = tsv.fields.select{|f| f.include?(" present") }
    removed = []
    tsv = tsv.select do |k,values|
      present = values.values_at(*present_fields).flatten.compact
      good = present.any? ||
        values["GO 981 terms"].any? || 
        (values["GO 43565 terms"].any? && (values["GO 6357 terms"].any?))

      removed << k unless good
      good
    end
    Open.write(file('removed'), removed * "\n")
    tsv
  end

  dep :filter_go
  task :process_go_fields => :tsv do
    tsv = step(:filter_go).load

    tsv.add_field "Has GO:0006355" do |k,values|
      values["GO 6355 terms"].any? ?  "Yes" : "No"
    end

    tsv.add_field "Has GO:0140223" do |k,values|
      values["GO 140223 terms"].any? ?  "Yes" : "No"
    end

    removed_go = ["GO 6355 terms", "GO 6355 evidence", "GO 140223 terms", "GO 140223 evidence"]
    good_fields = tsv.fields - removed_go

    tsv.slice good_fields
  end

  dep :process_go_fields
  task :cleanup => :tsv do
    tsv = dependencies.first.load

    tsv.select{|k,v| e = v["Entrez Gene ID"]; p = v["UniProt/SwissProt Accession"]; ! (e.nil? || e.empty?) &&  ! (p.nil? || p.empty?)}
  end

  task :orthologs => :tsv do
    organisms = %w(Hsa Mmu Rno)
    organisms.each do |target|
      organisms.each do |source|
        next if target == source
        Organism.send("ortholog_#{target}", source).produce
      end
    end
  end

end

#require 'TFCheckpoint/tasks/basic.rb'

#require 'rbbt/knowledge_base/TFCheckpoint'
#require 'rbbt/entity/TFCheckpoint'

