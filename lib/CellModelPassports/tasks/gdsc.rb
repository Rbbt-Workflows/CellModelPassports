
module CellModelPassports

  task :gdsc_AUC => :tsv do
    tsv = GDSC.drug_AUC.tsv
    tsv.key_field = "GDSC1000 drug ids"
    tsv = tsv.change_key "CCLE name", :identifiers => GDSC.drug_equivalences
    tsv = tsv.swap_id "COSMIC cell line ID", "GDSC1000 name", :identifiers => GDSC.cell_line_equivalences
    index = Association.index(tsv, :persist => false)
    index.add_field "Drug" do |k,v|
      k.split("~").first
    end
    index.add_field "Cell line" do |k,v|
      k.split("~").last
    end

    fields = index.fields[-2..-1] + index.fields[0..-3]
    index.reorder :key, fields
  end
end

