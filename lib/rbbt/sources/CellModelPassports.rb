require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/tsv/csv'

module CellModelPassports
  extend Resource
  self.subdir = 'share/databases/CellModelPassports'

  def self.organism(org="Hsa")
    "Hsa/feb2021"
  end

  def self.load_csv(url)
    TSV.open Open.open(url), :merge => true, :type => :double, :sep =>',', :header_hash => '', :zipped => true, :monitor => true, :namespace => CellModelPassports.organism
  end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  self.claim self.model_list, :csv,  "https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv"

  self.claim self.model_datasets, :csv,  "https://cog.sanger.ac.uk/cmp/download/model_dataset_availability_20220713.csv"

  self.claim self.driver_genes, :csv,  "https://cog.sanger.ac.uk/cmp/download/driver_genes_latest.csv"

  self.claim self.driver_gene_mutations, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/driver_mutations_latest.csv"
    tsv = TSV.csv Open.open(url), :merge => true, :type => :double, :namespace => CellModelPassports.organism
    enst2ensp = Organism.transcripts(CellModelPassports.organism).tsv :fields => ["Ensembl Protein ID"], :type => :single, :persist => true
    pos = ["ensembl_transcript", "ref_aa", "pep_coord", "alt_aa_list"].collect{|f| tsv.fields.index f }
    tsv.add_field "Mutated Isoform" do |k,values|
      Misc.zip_fields(values.values_at(*pos)).collect do |transcript,ref,pos,alt|
        next if alt.nil?
        protein = enst2ensp[transcript] || transcript
        alt.split(",").collect do |a|
          change = [ref, pos, a] * ""
          [protein, change] * ":"
        end.compact * ";"
      end
    end
  end

  self.claim self.predisposition_variants, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/cancer_predisposition_variants_latest.csv"
    tsv = TSV.csv Open.open(url), :merge => true, :type => :double, :namespace => CellModelPassports.organism
    enst2ensp = Organism.transcripts(CellModelPassports.organism).tsv :fields => ["Ensembl Protein ID"], :type => :single, :persist => true
    pos = ["ensembl_transcript", "ref_aa", "pep_coord", "alt_aa_list"].collect{|f| tsv.fields.index f }
    tsv.add_field "Mutated Isoform" do |k,values|
      Misc.zip_fields(values.values_at(*pos)).collect do |transcript,ref,pos,alt|
        next if alt.nil?
        protein = enst2ensp[transcript] || transcript
        alt.split(",").collect do |a|
          change = [ref, pos, a] * ""
          [protein, change] * ":"
        end.compact * ";"
      end
    end
  end

  self.claim self.driver_mutations, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/mutations_summary_latest.csv.gz"
    tsv = TSV.csv Open.open(url), :merge => true, :type => :double, :namespace => CellModelPassports.organism
  end

  self.claim self.all_mutations, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/mutations_all_latest.csv.gz"
    load_csv(url)
  end
   
  self.claim self.driver_cnv, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/cnv_summary_latest.csv.gz"
    TSV.csv Open.open(url), :merge => true, :type => :double, :namespace => CellModelPassports.organism
  end
   

  self.claim self.fusions, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/fusions_latest.csv.gz"
    load_csv(url)
  end

  self.claim self.proteomics, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/proteomics_latest.csv.gz"
    load_csv(url)
  end

  self.claim self.abundance, :proc do
    prot = CellModelPassports.proteomics.tsv(:key_field => 'symbol', :fields => ["model_id", "protein_intensity"])
    prot.key_field = "Associated Gene Name"
    models = prot.column('model_id').values.flatten.sort.uniq

    tsv = TSV.setup({}, :key_field => "Associated Gene Name", :fields => models, :cast => :to_f, :type => :list, :namespace => CellModelPassports.organism)

    TSV.traverse prot, :into => tsv, :bar => true do |gene,values|
      intensities = [nil] * models.length
      Misc.zip_fields(values).each do |model,intensity|
        intensities[models.index(model)] = intensity
      end
      [gene, intensities]
    end

    tsv.transpose("model_id")
  end



  self.claim self.gene_wes_cnv, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/WES_pureCN_CNV_genes_latest.csv.gz"
    load_csv(url)
  end

  self.claim self.expression, :proc do
    require 'rbbt/sources/organism'

    url = "https://cog.sanger.ac.uk/cmp/download/rnaseq_latest.csv.gz"

    models = {}
    symbols = Set.new
    TSV.traverse Open.open(url), :type => :list, :bar => true do |line| 
      next if line =~ /^dataset_id/
      dataset_id, id, model_id, gene_id, read_count, fpkm, tpm, data_source, duplicate, dataset_name, model_name, symbol = parts = line.split(",")
      model_values = models[model_id] ||= {}
      symbols << symbol
      model_values[symbol] ||= []
      model_values[symbol] << [fpkm, tpm,read_count,data_source]
    end

    symbols = symbols.to_a.sort

    fpkm = TSV.setup({}, :key_field => "model_id", :fields => symbols, :type => :list, :cast => :to_f, :namespace => CellModelPassports.organism)

    models.each do |model_id,values|
      symbol_values = values.chunked_values_at(symbols)

      symbol_tpm = symbol_values.collect do |l| 
        if l.nil?
          []
        else
          l.collect{|a| a[1].to_f }
        end
      end

      fpkm[model_id] = symbol_tpm
    end
  
    fpkm
  end
end

if __FILE__ == $0
  Log.severity = 0
  methods = CellModelPassports.resources.keys.collect{|p| p.split("/").last }

  #methods.each do |method|
  #  CellModelPassports[method].produce
  #end

  file = CellModelPassports[methods.last].produce
  ppp CellModelPassports[methods.last].tsv.head.to_s
  Log.tsv CellModelPassports[methods.last].tsv
  
end

