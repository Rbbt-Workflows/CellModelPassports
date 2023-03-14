module CellModelPassports
  task :expression_for_decoupler => :tsv do
    tsv = CellModelPassports.expression.tsv
    tsv.each do |k,vs|
      tsv[k] = vs.collect{|v| v.nil? ? 0 : v }
    end
    tsv
  end

  dep SaezLab, :regulome, :jobname => 'Default'
  dep :expression_for_decoupler
  dep_task :tf_activity, SaezLab, :decoupler, :matrix => :expression_for_decoupler, :network =>  :regulome


  dep :cell_model_expression
  task :expression_for_decoupler_by_cell_model => :tsv do
    cell_model = step(:cell_model_expression).recursive_inputs[:cell_model]
    tsv = step(:cell_model_expression).load.transpose("model_id")
    tsv.cast = :to_f
    tsv[cell_model] = tsv.delete(tsv.keys.first)
    tsv
  end

  dep SaezLab, :regulome, :jobname => 'Default'
  dep :expression_for_decoupler_by_cell_model
  dep_task :cell_model_decoupler, SaezLab, :decoupler, :matrix => :expression_for_decoupler_by_cell_model, :network =>  :regulome, :_method => :consensus

  dep :cell_model_decoupler
  input :threshold, :float, "Decoupler p-value threshold", 0.1
  task :cell_model_tf_activity => :tsv do |threshold|
    tsv = step(:cell_model_decoupler).load
    tsv = tsv.transpose("Associated Gene Name").to_single
    genes = tsv.keys
    p_values = genes.select{|g| g.include? 'pvalue'}
    genes -= p_values
    p_values.each do |field|
      p = tsv.delete field
      gene = field.sub(" (pvalue)", '')
      tsv.delete gene if p.to_f > threshold
    end
    tsv.fields = ["TF Activity"]
    tsv
  end

end
