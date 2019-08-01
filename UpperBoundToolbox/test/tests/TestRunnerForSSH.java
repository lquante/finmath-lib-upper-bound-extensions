package tests;

import net.finmath.exception.CalculationException;

public class TestRunnerForSSH {
	
	// main parameters to adjust for simulation
	static int numberOfPaths = 100000;
	static int numberOfSubsimulationPathsA = 10000;
	static int numberOfSubsimulationPathsB = 10000;
	
	static int numberOfPeriods = 10;
	public  static void main(String[] args) throws CalculationException {
		// set parameters
		CreateTestModel.setNumberOfPaths(numberOfPaths);
		CreateTestBermudanSwaption.setNumberOfExercisePeriods(numberOfPeriods);
		TestValuationMethods.setNumberOfSubsimulationsStepA(numberOfSubsimulationPathsA);
		TestValuationMethods.setNumberOfSubsimulationsStepB(numberOfSubsimulationPathsB);
	
		// run test
		
		TestValuationMethods.testSwaptionValuationMethod();
	
	}

}