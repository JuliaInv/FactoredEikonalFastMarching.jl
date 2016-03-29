using Base.Test

allPassed = true
try 
	include("sensitivityTest.jl");
catch
	allPassed = false
	warn("including sensitivityTest.jl failed.")
end

try 
	runSensitivityTest();
catch
	allPassed = false
	warn("runSensitivityTest() had test errors")
end

@test allPassed == true
