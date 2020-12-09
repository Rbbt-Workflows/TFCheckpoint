require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/tsv/excel'
require 'rbbt/sources/organism'

module TFCheckpoint
  extend Resource
  self.subdir = 'share/databases/TFCheckpoint2'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib



  TFCheckpoint.claim TFCheckpoint.Saere, :proc do 
    tsv = TSV.xlsx Rbbt.data["Saere_et_al.xls"], :skip_rows => 1, :merge => true, :key_field => 1
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Hsa"
    tsv.key_field = "Associated Gene Name"
    tsv.to_single
  end

  TFCheckpoint.claim TFCheckpoint.Ravazi, :proc do 
    tsv = TSV.excel Rbbt.data["Ravazi.xls"]
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Hsa"
    tsv.key_field = "Associated Gene Name"
    tsv = tsv.select{|k,v| ! k.include? "These listings"}
    tsv.to_single
  end

  TFCheckpoint.claim TFCheckpoint.TcoF_human, :proc do 
    tsv = TSV.excel Rbbt.data["TcoF-DB_human.xlsx"], :key_field => 1, :fields => []
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Hsa"
    tsv.key_field = "Associated Gene Name"
    tsv.to_single
  end

  TFCheckpoint.claim TFCheckpoint.TcoF_mouse, :proc do 
    tsv = TSV.excel Rbbt.data["TcoF-DB_mouse.xlsx"], :key_field => 1, :fields => []
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Mmu"
    tsv.key_field = "Associated Gene Name"
    tsv.to_single
  end

  TFCheckpoint.claim TFCheckpoint.animal_tfdb_Mus_musculus, :proc do 
    tsv = TSV.open Rbbt.data["animal_tfdb_Mus_musculus.txt"], :header_hash => '', :key_field => "Symbol", :fields => [], :type => :single
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Mmu"
    tsv.key_field = "Associated Gene Name"
    tsv.to_single
  end

  TFCheckpoint.claim TFCheckpoint.animal_tfdb_Rattus_norvegicus, :proc do 
    tsv = TSV.open Rbbt.data["animal_tfdb_Rattus_norvegicus.txt"], :header_hash => '', :key_field => "Symbol", :fields => [], :type => :single
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Rno"
    tsv.key_field = "Associated Gene Name"
    tsv.to_single
  end

  TFCheckpoint.claim TFCheckpoint.animal_tfdb_homo_sapiens, :proc do 
    tsv = TSV.open Rbbt.data["animal_tfdb_homo_sapiens.txt"], :header_hash => '', :key_field => "Symbol", :fields => [], :type => :single
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Hsa"
    tsv.key_field = "Associated Gene Name"
    tsv.to_single
  end

  %w(go_43565_human go_43565_mouse go_43565_rat go_6357_human go_6357_mouse
  go_6357_rat go_981_human go_981_mouse go_981_rat).each do |db|
    namespace = case db.split("_").last
                when 'human'
                  'Hsa'
                when 'mouse'
                  'Mmu'
                when 'rat'
                  'Rno'
                end
    TFCheckpoint.claim TFCheckpoint[db], :proc do 
      tsv = TSV.open Rbbt.data["#{db}.tsv"], :header_hash => '', :key_field => 1, :fields => [4,7], :type => :double, :merge => true
      tsv = tsv.reorder :key, [:key] + tsv.fields, :one2one => true, :merge => true
      tsv.namespace = namespace
      tsv.key_field = "UniProt/SwissProt Accession"
      tsv = tsv.change_key "Associated Gene Name"
      tsv
    end
  end

  TFCheckpoint.claim TFCheckpoint.lambert_2018, :proc do 
    tsv = TSV.open Rbbt.data["lambert_2018.txt"], :header_hash => '', :key_field => "HGNC symbol", :fields => ["Is TF?"], :type => :list
    tsv = tsv.reorder :key, [:key, "Is TF?"]
    tsv.namespace = "Hsa"
    tsv.key_field = "Associated Gene Name"
    tsv
  end

  TFCheckpoint.claim TFCheckpoint.vaquerizas, :proc do 
    tsv = TSV.open Rbbt.data["vaquerizas.txt"], :header_hash => '', :key_field => 5, :fields => [1,0], :type => :list
    tsv = tsv.reorder :key, [:key, "Class"]
    tsv.namespace = "Hsa"
    tsv.key_field = "Associated Gene Name"
    tsv
  end

  ["human", "mouse", "rat"].each do |org|
    namespace = case org
                when 'human'
                  'Hsa'
                when 'mouse'
                  'Mmu'
                when 'rat'
                  'Rno'
                end
    TFCheckpoint.claim TFCheckpoint["jaspar_" + org], :proc do 
      list = Rbbt.data["jaspar_#{org}.txt"].list.collect{|e| e.split(";")}.flatten.collect{|e| e.strip}.reject{|e| e.empty?}
      tsv = TSV.setup(list, :key_field => "UniProt/SwissProt Accession", :fields => [], :type => :list, :namespace => namespace)
      tsv.attach Organism.identifiers(namespace), :fields => ["Associated Gene Name"]
      tsv = tsv.select("Associated Gene Name"){|n| ! n.nil? && ! n.empty?}
      tsv = tsv.reorder "Associated Gene Name", [:key]
      tsv
    end
  end

  # UPDATE from Oct2020
  
  TFCheckpoint.claim TFCheckpoint.humantfdb, :proc do 
    tsv = TSV.open Rbbt.data["humantfdb.txt"], :header_hash => '', :key_field => 1, :fields => [3], :type => :single
    tsv = tsv.reorder :key, [:key]
    tsv.namespace = "Hsa"
    tsv.key_field = "Associated Gene Name"
    tsv
  end

  TFCheckpoint.claim TFCheckpoint.GREEKC, :proc do 
    tsv = TSV.setup(Rbbt.data["dbTF_gene_product_set.tsv"].read.split("\n").collect{|l| l.split("\t").first}, "Associated Gene Name~#:type=:single")
    tsv.namespace = "Hsa"
    tsv = tsv.reorder :key, [:key]
    tsv
  end

  TFCheckpoint.claim TFCheckpoint.DBD, :proc do 
    tsv = TSV.open Rbbt.data["dbd_orfeome_tfcat.tsv"], :header_hash => '', :fields => [3], :type => :double
    tsv = tsv.select{|k,v| v != "0"}
    tsv.namespace = "Hsa"
    tsv.key_field = "Entrez Gene ID"
    tsv = tsv.reorder :key, [:key]
    tsv
  end

  TFCheckpoint.claim TFCheckpoint.ORFeome, :proc do 
    tsv = TSV.open Rbbt.data["dbd_orfeome_tfcat.tsv"], :header_hash => '', :fields => [4], :type => :list
    tsv = tsv.select{|k,v| v != "0"}
    tsv.key_field = "Entrez Gene ID"
    tsv.namespace = "Hsa"
    tsv = tsv.reorder :key, [:key]
    tsv
  end

  TFCheckpoint.claim TFCheckpoint.TFCat, :proc do 
    tsv = TSV.open Rbbt.data["dbd_orfeome_tfcat.tsv"], :header_hash => '', :fields => [5], :type => :list
    tsv = tsv.select{|k,v| v != "0"}
    tsv.key_field = "Entrez Gene ID"
    tsv.namespace = "Hsa"
    tsv = tsv.reorder :key, [:key]
    tsv
  end
end

if __FILE__ == $0
  update = true

  Log.tsv TFCheckpoint.humantfdb.produce(true).tsv
  Log.tsv TFCheckpoint.GREEKC.produce(true).tsv
  Log.tsv TFCheckpoint.DBD.produce(true).tsv
  Log.tsv TFCheckpoint.ORFeome.produce(true).tsv
  Log.tsv TFCheckpoint.TFCat.produce(true).tsv

  #if update

  #  TFCheckpoint.resources.keys.each do |key|
  #    TFCheckpoint[key.split("/").last].produce(false).tsv
  #  end 

  #else

  #  file = TFCheckpoint.resources.keys.last.split("/").last
  #  Log.tsv TFCheckpoint[file].produce(true).tsv
  #end
end

