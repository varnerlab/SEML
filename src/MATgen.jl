using MAT
# using PyCall

function generate_MAT_file(rnx::Array, species::Array, stoichiometry::Array)
  filename = "model.mat"
  file = matopen(filename, "w")
  rxns = [rr.rnxName for rr in rnx]
  write(file, "rxns", rxns)
  write(file, "S", stoichiometry)
  lb = zeros(length(rxns), 1)
  write(file, "lb", lb)
  ub = 1000 * ones(length(rxns), 1)
  write(file, "ub", ub)
  c = zeros(length(rxns), 1)
  write(file, "c", c)
  write(file, "mets", species)
  b = 100 * ones(length(species), 1)
  write(file, "b", b)
  write(file, "osenseStr", "max")
  csense = ["L" for i = 1:length(species)]
  write(file, "csense", csense)

  # write(file, "genes", [""])
  # write(file, "rules", ["" for i = 1:length(rxns)])
  close(file)
  return filename
end
