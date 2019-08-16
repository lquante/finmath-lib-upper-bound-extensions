package tests;

import org.junit.Test;

import net.finmath.exception.CalculationException;

public class TestRunnerForSSH {
	
	// main parameters to adjust for simulation
	static int numberOfPaths = 10000;
	static int numberOfSubsimulationPathsA = 1000;
	static int numberOfSubsimulationPathsB = 1000;
	
	static int numberOfPeriods = 10;
	@Test
	public void testWithManyPaths() throws CalculationException
	{
		main(null);
	}
	
	public  static void main(String[] args) throws CalculationException {
		// set parameters
		CreateTestModel.setNumberOfPaths(numberOfPaths);
		CreateTestBermudanSwaption.setNumberOfExercisePeriods(numberOfPeriods);
		TestValuationMethods.setNumberOfSubsimulationsStepA(numberOfSubsimulationPathsA);
		TestValuationMethods.setNumberOfSubsimulationsStepB(numberOfSubsimulationPathsB);
	
		// run test
		
		TestValuationMethods.executePrintSwaptionValuationMethods();
	
	}

}
