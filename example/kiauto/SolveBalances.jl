include("./Include.jl")

function SolveBalances(TStart, TStop, TStep, dataDictionary)
	t = collect(TStart:TStep:TStop)
	all_species_reversed_dict = dataDictionary["all_species_reversed_dict"]
	# initial species concentration
	y0 = dataDictionary["initial_condition"]
	# call ODE solver
	f(t, y) = Balances(t, y, dataDictionary)
	t, y = ode23s(f, y0, t; points=:specified)

	# data transfer
	row = length(t)
	col = length(y0)
	Y = zeros(row, col)
	foreach(x->(Y[x, :] = y[x]), collect(1:row))
	# plot results
	plt.figure("simulation $TStart to $TStop")
	subplt_col = ceil(Int, col/4)
	for i = 1:col
		plt.subplot(4, subplt_col, i)
		plt.plot(t, Y[:,i])
		plt.title(all_species_reversed_dict[i])
	end
	plt.show()

	return (t, Y)
end


# # uncomment to test
# dataDictionary = generate_model_parameters_dictionary()
# SolveBalances(0, 10, 0.1, dataDictionary)
