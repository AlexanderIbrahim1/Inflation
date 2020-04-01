module Inflation

export convert_01tuples
export rows_and_cols
export get_M_matrix
export get_b_array
export get_v_array_Variable
export integrate_pfunc

export pfunc1
export pfunc2
export get_pfunct

using Convex
using GLPK
using Printf

export InflKind
export cut
export spiral

@enum InflKind cut spiral

function convert_01tuples(str::String, char::Char)::Vector{String}
	# accepts a string with a certain character 'char'
	# returns all possible combinations where 'char' is
	# replaced with either 0 or 1
	
	foundall  = false
	newcombos = Vector{String}([str])
	combos    = Vector{String}([])
	while !foundall
		foundall  = true
		combos    = newcombos
		newcombos = []
		for combo in combos
			pos = findfirst(char, combo)
			if pos != nothing
				combo0 = replace(combo, char=>'0', count=1)
				combo1 = replace(combo, char=>'1', count=1)
				push!(newcombos, combo0)
				push!(newcombos, combo1)
				foundall = false
			end
		end
	end
	
	return combos
end

function rows_and_cols(infl_strings::Vector{String})::Tuple{Int64, Int64}
	# returns the number of rows and columns of the M-vector
	if isempty(infl_strings)
		error("'infl_strings' must have at least one element")
	end
	
	str1 = infl_strings[1]
	infl_strings_y = convert_01tuples(str1, 'y')
	
	nrows = length(infl_strings)*length(infl_strings_y)
	ncols = 2^length(str1)
	
	return (nrows, ncols)
end

function get_M_matrix(infl_strings::Vector{String})::Array{Int64, 2}
	# returns the M-vector needed for the inflation algorithm
	# infl_strings : array of strings that determine the inflation combinations
	# EXAMPLE: infl_strings = ['yyx', 'yxy', 'xyy']
	
	Nrows, Ncols = rows_and_cols(infl_strings)
	Mvec = zeros(Int64, Nrows, Ncols)
	
	for (it, str_yx) in enumerate(infl_strings)
		infl_strings_x = convert_01tuples(str_yx, 'y')
		N_x            = length(infl_strings_x)
		for (ix, str_x) in enumerate(infl_strings_x)
			i = (ix-1) + N_x*(it-1)
			fill_row_Mvec(i, str_x, Mvec)
		end
	end
	
	return Mvec
end

function fill_row_Mvec(i::Int64, str_x::String, Mvec::Array{Int64, 2})::Nothing
	# fills the i^th row of Mvec using the combinations given by str_x
	for str in convert_01tuples(str_x, 'x')
		ic = parse(Int64, str, base=2)
		Mvec[i+1, ic+1] = 1
	end
end

function get_b_array(infl_strings::Vector{String}, pfunc::Function, kind::InflKind)::Array{Float64, 1}
	# returns the b-vector for the Cut Inflation of the triangle causal model
	[@assert only_contains("xy", str) for str in infl_strings]
	
	Nrows_per_key = length(convert_01tuples(infl_strings[1], 'y'))
	Nkeys         = length(infl_strings)
	barr          = Vector{Float64}(undef, Nrows_per_key*Nkeys)
	
	for (i, str_yx) in enumerate(infl_strings)
		infl_strings_x = convert_01tuples(str_yx, 'y')
		for (j, str_x) in enumerate(infl_strings_x)
			#index = (i-1) + Nkeys*(j-1) + 1
			index = (j-1) + Nrows_per_key*(i-1) + 1
			barr[index] = inflation_RHS(str_x, pfunc, kind)
		end
	end
	
	return barr
end

function get_v_array_Variable(infl_strings::Vector{String})
	# gives the Variable instance with the correct size for the v-array
	Nv_size = 2^length(infl_strings[1])
	return Variable(Nv_size)
end

function inflation_RHS_Cut(str_x::String, pfunc::Function)::Float64
	# performs the RHS integration for the Cut inflation
	
	# check if 'str_x' is valid for the Cut inflation integration
	Nx = sum([char == 'x' for char in str_x])
	@assert Nx            == 1
	@assert length(str_x) == 3
	
	total = 0.0
	if str_x[3] == 'x'
		# eg) '00x', or '01x', ...
		strA = str_x[1] * "xx"
		strB = "x" * str_x[2] * "x"
		totalA = integrate_pfunc(pfunc, strA, 'x')
		totalB = integrate_pfunc(pfunc, strB, 'x')
		total += totalA*totalB
	else
		# eg) 'x00', '0x1', ...
		total += integrate_pfunc(pfunc, str_x, 'x')
	end
	
	return total
end

function inflation_RHS_Spiral(str_x::String, pfunc::Function)::Float64
	# performs the RHS integration for the Spiral inflation
	
	# check if 'str_x' is valid for the Spiral inflation integration
	Nx = sum([char == 'x' for char in str_x])
	@assert Nx            == 3
	@assert length(str_x) == 6
	
	total = 0.0
	if     match_str(str_x, [2, 4, 6], 'x')
		# 'yxyxyx'
		str = str_x[1] * str_x[3] * str_x[5]
		total += pfunc(str)
	elseif match_str(str_x, [2, 3, 5], 'x')
		# 'yxxyxy'
		str   = str_x[1] * str_x[4] * str_x[6]
		strAB = str[1:2] * "x"
		strC  = "xx" * str[3]
		
		totalAB = integrate_pfunc(pfunc, strAB, 'x')
		totalC  = integrate_pfunc(pfunc, strC,  'x')
		total  += totalAB*totalC
	elseif match_str(str_x, [1, 4, 5], 'x')
		# 'xyyxxy'
		str   = str_x[2] * str_x[3] * str_x[6]
		strBC = "x" * str[2:3]
		strA  = str[1] * "xx"
		
		totalBC = integrate_pfunc(pfunc, strBC, 'x')
		totalA  = integrate_pfunc(pfunc, strA,  'x')
		total  += totalA*totalBC
	elseif match_str(str_x, [1, 3, 6], 'x')
		# 'xyxyyx'
		str   = str_x[2] * str_x[4] * str_x[5]
		strAC = str[1] * "x" * str[3]
		strB  = "x" * str[2] * "x"
		
		totalAC = integrate_pfunc(pfunc, strAC, 'x')
		totalB  = integrate_pfunc(pfunc, strB,  'x')
		total  += totalAC*totalB
		
	elseif match_str(str_x, [1, 3, 5], 'x')
		# 'xyxyxy'
		str  = str_x[2] * str_x[4] * str_x[6]
		strA = str[1] * "xx"
		strB = "x" * str[2] * "x"
		strC = "xx" * str[3]
		
		totalA = integrate_pfunc(pfunc, strA, 'x')
		totalB = integrate_pfunc(pfunc, strB, 'x')
		totalC = integrate_pfunc(pfunc, strC, 'x')
		total += totalA*totalB*totalC
	else
		error("Invalid inflation string.")
	end
	
	return total
end

function match_str(str_x::String, indices::Vector{Int64}, char::Char)::Bool
	# check if 'str_x' has 'char' at certain indices
	matches = true
	for i in indices
		if str_x[i] != char
			matches = false
			break
		end
	end
	
	return matches
end

function inflation_RHS(str_x::String, pfunc::Function, kind::InflKind)::Float64
	# select which kind of RHS inflation to use
	if     kind == cut
		return inflation_RHS_Cut(str_x, pfunc)
	elseif kind == spiral
		return inflation_RHS_Spiral(str_x, pfunc)
	else
		error("Only the inflations for 'Cut' and 'Spiral' have been implemented.")
	end
end

function integrate_pfunc(pfunc::Function, str_x::String, char::Char)::Float64
	# pfunc : the probability distribution for the triangle causal model
	# str_x : a string of 0's, 1's, and 'char's
	#       : the 'char's will be replaced with 0's and 1's
	#       : and the combination of all these terms will be added together
	
	# check if str_x only has 0's, 1's, and char's
	@assert only_contains("01"*char, str_x)
	
	infl_strings_01 = convert_01tuples(str_x, char)
	total = sum([pfunc(str) for str in infl_strings_01])
	
	return total
end

function get_pfunct(t::Float64)::Function
	# for problem 2, return pfunc(a, b, c; t)
	@assert 0.0 <= t <= 1.0
	
	function pfunct(str::String)::Float64
		# returns the 1st probability distribution for the triangle causal model
		@assert (a = tryparse(Int64, string(str[1]))) != nothing
		@assert (b = tryparse(Int64, string(str[2]))) != nothing
		@assert (c = tryparse(Int64, string(str[3]))) != nothing
		@assert (a==0 || a==1)
		@assert (b==0 || b==1)
		@assert (c==0 || c==1)
		
		if (a == b == c)
			return t*0.5
		else
			return (1.0 - t)/6.0
		end
	end
	
	return pfunct
end

function pfunc1(str::String)::Float64
	@assert length(str) == 3
	@assert (a = tryparse(Int64, string(str[1]))) != nothing
	@assert (b = tryparse(Int64, string(str[2]))) != nothing
	@assert (c = tryparse(Int64, string(str[3]))) != nothing
	return pfunc1(a, b, c)
end

function pfunc1(a::Int64, b::Int64, c::Int64)::Float64
	# returns the 1st probability distribution for the triangle causal model
	@assert (a==0 || a==1)
	@assert (b==0 || b==1)
	@assert (c==0 || c==1)
	
	if (a == b == c)
		return 0.5
	else
		return 0.0
	end
end

function pfunc2(str::String)::Float64
	@assert length(str) == 3
	@assert (a = tryparse(Int64, string(str[1]))) != nothing
	@assert (b = tryparse(Int64, string(str[2]))) != nothing
	@assert (c = tryparse(Int64, string(str[3]))) != nothing
	return pfunc2(a, b, c)
end

function pfunc2(a::Int64, b::Int64, c::Int64)::Float64
	# returns the 2nd probability distribution for the triangle causal model
	@assert (a==0 || a==1)
	@assert (b==0 || b==1)
	@assert (c==0 || c==1)
	
	if (a + b + c == 1)
		return 1.0/3.0
	else
		return 0.0
	end
end

function only_contains(ok_chars::String, str::String)::Bool
	# Checks that 'str' only contains elements from 'ok_chars'
	return all([c in ok_chars for c in str])
end

end # module
