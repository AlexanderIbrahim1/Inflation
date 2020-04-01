using Inflation
using Printf
using Test

@testset "convert_01tuples" begin
	# Test 1
	str  = "0xx"
	char = 'x'
	exp_output = ["000", "001", "010", "011"]
	act_output = convert_01tuples(str, char)
	@test exp_output == act_output
	
	# Test 2
	str  = "01x"
	char = 'x'
	exp_output = ["010", "011"]
	act_output = convert_01tuples(str, char)
	@test exp_output == act_output
	
	# Test 3
	str  = "y1x0"
	char = 'y'
	exp_output = ["01x0", "11x0"]
	act_output = convert_01tuples(str, char)
	@test exp_output == act_output
end;

@testset "rows_and_cols" begin
	# Test 1
	infl_strings = ["yyx", "yxy", "xyy"]
	(a, b) = rows_and_cols(infl_strings)
	@test (a, b) == (12, 8)
end;

@testset "get_M_matrix" begin
	# Test 1
	infl_strings = ["yyx", "yxy", "xyy"]
	exp_output   = [1 1 0 0 0 0 0 0;
					0 0 1 1 0 0 0 0;
					0 0 0 0 1 1 0 0;
					0 0 0 0 0 0 1 1;
					1 0 1 0 0 0 0 0;
					0 1 0 1 0 0 0 0;
					0 0 0 0 1 0 1 0;
					0 0 0 0 0 1 0 1;
					1 0 0 0 1 0 0 0;
					0 1 0 0 0 1 0 0;
					0 0 1 0 0 0 1 0;
					0 0 0 1 0 0 0 1]
	act_output = get_M_matrix(infl_strings)
	@test exp_output == act_output
end;

@testset "pfunc_series" begin
	# Testing pfunc1
	@test pfunc1(0, 0, 0) == 0.5
	@test pfunc1(0, 0, 1) == 0.0
	@test pfunc1(0, 1, 0) == 0.0
	@test pfunc1(0, 1, 1) == 0.0
	@test pfunc1(1, 0, 0) == 0.0
	@test pfunc1(1, 0, 1) == 0.0
	@test pfunc1(1, 1, 0) == 0.0
	@test pfunc1(1, 1, 1) == 0.5
	
	@test pfunc1("000") == 0.5
	@test pfunc1("001") == 0.0
	
	# Testing pfunc2
	@test pfunc2(0, 0, 0) == 0.0
	@test pfunc2(0, 0, 1) == 1.0/3.0
	@test pfunc2(0, 1, 0) == 1.0/3.0
	@test pfunc2(0, 1, 1) == 0.0
	@test pfunc2(1, 0, 0) == 1.0/3.0
	@test pfunc2(1, 0, 1) == 0.0
	@test pfunc2(1, 1, 0) == 0.0
	@test pfunc2(1, 1, 1) == 0.0
	
	@test pfunc2("000") == 0.0
	@test pfunc2("001") == 1.0/3.0
	
	# Testing pfunct
	pfunct = get_pfunct(0.5)
	@test pfunct("000") == 1.0/4.0
	@test pfunct("100") == 1.0/12.0
end;

@testset "integrate_pfunc" begin
	# Test 1 : pfunc1(0, 0, 0) + pfunc1(0, 0, 1) == 0.5 + 0.0 == 0.5
	@test integrate_pfunc(pfunc1, "00x", 'x') == 0.5
	
	# Test 2 : pfunc1(0, 1, 1) + pfunc1(0, 0, 1) == 0.0 + 0.0 == 0.0
	@test integrate_pfunc(pfunc1, "0x1", 'x') == 0.0
	
	# Test 3 :   pfunc1(0, 1, 1) + pfunc1(0, 0, 1) 
	#        : + pfunc1(0, 0, 0) + pfunc1(0, 1, 0)
	#        : = 0.0 + 0.0 + 0.5 + 0.0 == 0.5
	@test integrate_pfunc(pfunc1, "0xx", 'x') == 0.5
end;

@testset "spiral_inflation_test" begin
	# The spiral inflation for the triangle causal model has a very
	# large M-matrix and b-array, so we can test a few terms
	
	infl_strings = ["yxyxyx", "yxxyxy", "xyyxxy", "xyxyyx", "xyxyxy"]
	
	# creating the M-matrix and b-array to be used in the Convex.jl optimization
	M = get_M_matrix(infl_strings)
	b = get_b_array(infl_strings, pfunc2, spiral)
	
	############################################################################
	# P_{A1,B1,C1}(0, 0, 0) term?
	############################################################################
	# What is the index for this term?
	# index :   8 * [(position in infl_strings) - 1] + [binary('000') + 1]
	# index : = 8*(0) + 1 = 1
	
	# What does the row in the M-matrix look like?
	# infl_string = '0x0x0x'
	# length = 2^6 = 64
	# '000000' => 1      '000001' => 2      '000100' => 5      '000101' => 6 
	# '010000' => 17     '010001' => 18     '010100' => 21     '010101' => 22
	Mvec000 = zeros(Int64, 64)
	Mvec000[1]  = 1  ;  Mvec000[2]  = 1  ;  Mvec000[5]  = 1  ;  Mvec000[6]  = 1
	Mvec000[17] = 1  ;  Mvec000[18] = 1  ;  Mvec000[21] = 1  ;  Mvec000[22] = 1
	
	# What does the b-array element look like?
	# P_{A1,B1,C1}(0,0,0) = P_{A,B,C}(0,0,0) = 0.0
	b000 = 0.0
	
	# check if correct
	@test M[1,:] == Mvec000
	@test b[1] == b000
	
	############################################################################
	# P_{A2,B2,C1}(1, 0, 0) term?
	############################################################################
	# What is the index for this term?
	# index :   8 * [(position in infl_strings) - 1] + [binary('100') + 1]
	# index : = 8*(3) + 5 = 29
	
	# What does the row in the M-matrix look like?
	# infl_string = 'x1x00x'
	# length = 2^6 = 64
	# '010000' => 17     '010001' => 18     '011000' => 25     '011001' => 26
	# '110000' => 49     '110001' => 50     '111000' => 57     '111001' => 58 
	Mvec100 = zeros(Int64, 64)
	Mvec100[17] = 1  ;  Mvec100[18] = 1  ;  Mvec100[25] = 1  ;  Mvec100[26] = 1
	Mvec100[49] = 1  ;  Mvec100[50] = 1  ;  Mvec100[57] = 1  ;  Mvec100[58] = 1
	
	# What does the b-array element look like?
	# P_{A1,B1,C1}(1,0,0) = P_{A,C}(1,0) P_{B}(0) = [0 + 1/3][0 + 1/3 + 1/3 + 0]
	#                     = 2/9 = 0.2222...
	b100 = 2.0/9.0
	
	@test M[29,:] == Mvec100
	@test b[29] == b100
end;

