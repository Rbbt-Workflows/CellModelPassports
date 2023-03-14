require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/matrix'
require 'rbbt/matrix/barcode'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/CellModelPassports'

Workflow.require_workflow "SaezLab"
Workflow.require_workflow "GDSC"

module CellModelPassports
  extend Workflow

  task :cell_models => :array do
    tsv = CellModelPassports.model_list.tsv :persist => true
    tsv.column("model_name").values.sort
  end

  helper :cell_model_id do |cell_model|
    index = CellModelPassports.model_list.index :persist => true
    index[cell_model]
  end

  helper :load_tsv do |file, fields|
    file.tsv :key_field => 'model_id', :fields => fields, :type => :double, :merge => true, :zipped => true, :monitor => true, :persist => :preload
  end

  helper :extract_model_tsv do |tsv,model_id|
    tsv = tsv.subset([model_id])
    if tsv.fields.include?("symbol")
      tsv = tsv.reorder "symbol", tsv.fields - ["symbol"], :zipped => true
    elsif tsv.fields.include?("gene_symbol")
      tsv = tsv.reorder "gene_symbol", tsv.fields - ["gene_symbol"], :zipped => true
    end
    tsv
  end


  input :cell_model, :string, "Cell model name", nil, :jobname => true
  task :cell_model_mutations => :tsv do |cell_model|

    model_id = cell_model_id cell_model

    tsv = load_tsv CellModelPassports.all_mutations, %w(gene_symbol protein_mutation)
    extract_model_tsv tsv, model_id
  end

  input :cell_model, :string, "Cell model name", nil, :jobname => true
  task :cell_model_cnv => :tsv do |cell_model|

    model_id = cell_model_id cell_model

    tsv = load_tsv CellModelPassports.gene_wes_cnv, %w(symbol total_copy_number minor_copy_number loh cn_category comment)
    extract_model_tsv tsv, model_id
  end

  input :cell_model, :string, "Cell model name", nil, :jobname => true
  task :cell_model_expression => :tsv do |cell_model|

    model_id = cell_model_id cell_model

    tsv = CellModelPassports.expression.tsv :persist => true, :monitor => true
    tsv = tsv.subset([model_id])
    tsv = tsv.transpose("Associated Gene Name").to_single
    tsv.each{|k,v| tsv.delete k if v.nil? }
    tsv.fields = ["TPM Expression"]
    tsv
  end

  helper :matrix_activity do |tsv_file|
    m = RbbtMatrix.new tsv_file
    m = m.transpose
    tsv = TSV.open(m.to_average.to_activity.data_file, :type => :list, :cast => :to_i)
    tsv.each do |k,v|
      fixed = v.collect do |e| 
        case e.to_i
        when 1
          -1
        when 2
          1
        else
          nil
        end
      end
      v.replace(fixed)
    end
    tsv
  end

  input :cell_model, :string, "Cell model name", nil, :jobname => true
  task :cell_model_expression_cluster => :tsv do |cell_model|

    model_id = cell_model_id cell_model

    begin
      tsv = matrix_activity CellModelPassports.expression
      raise ParameterException, "#{ model_id } field not found in matrix" unless tsv.fields.include?(model_id)
      tsv = tsv.slice(model_id).to_single
      tsv.delete_if{|k,v| v.nil? }
      tsv.fields = ["Transcription Level Cluster"]
      tsv
    rescue
      Log.exception $!
      TSV.setup({}, "Associated Gene Name~Transcription Level Cluster#:type=:single#:cast=:to_i")
    end
  end


  input :cell_model, :string, "Cell model name", nil, :jobname => true
  task :cell_model_proteomics => :tsv do |cell_model|

    model_id = cell_model_id cell_model

    tsv = load_tsv CellModelPassports.proteomics, %w(symbol uniprot_id protein_intensity zscore)
    extract_model_tsv tsv, model_id
  end

  input :cell_model, :string, "Cell model name", nil, :jobname => true
  task :cell_model_proteomics_cluster => :tsv do |cell_model|

    model_id = cell_model_id cell_model

    begin
      tsv = matrix_activity CellModelPassports.abundance
      raise ParameterException, "#{ model_id } field not found in matrix" unless tsv.fields.include?(model_id)
      tsv = tsv.slice(model_id).to_single
      tsv.delete_if{|k,v| v.nil? }
      tsv.fields = ["Protein Intensity Level Cluster"]
      tsv
    rescue
      Log.exception $!
      TSV.setup({}, "Associated Gene Name~Protein Intensity Level Cluster#:type=:single#:cast=:to_i")
    end
  end


  dep :cell_model_mutations
  dep :cell_model_cnv
  dep :cell_model_expression
  dep :cell_model_expression_cluster
  dep :cell_model_tf_activity
  dep :cell_model_proteomics
  dep :cell_model_proteomics_cluster
  task :cell_model_bundle => :tsv do
    dependencies.inject(nil) do |acc,d|
      next acc if d.error?
      tsv = d.load.to_double
      tsv.key_field = "Associated Gene Name"
      if acc.nil?
        acc = tsv
      else
        acc.attach tsv, :complete => true
      end
      acc
    end
  end

end

require 'CellModelPassports/tasks/decopler.rb'
require 'CellModelPassports/tasks/gdsc.rb'

#require 'rbbt/knowledge_base/CellModelPassports'
#require 'rbbt/entity/CellModelPassports'

