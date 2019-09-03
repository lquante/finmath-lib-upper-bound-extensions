package tests;

import org.junit.Test;

import net.finmath.exception.CalculationException;

public class TestRunnerForSSH {
	// model parameters
	static double lastTimePoint = 30;
	static double timeDiscretizationLength = 0.25;
	static double liborPeriodLength = 0.5;
	// monte carlo parameters
	static int numberOfPaths = 10000;
	static int numberOfSubsimulationsStepA = 10;
	static int numberOfSubsimulationsStepB = 10;
	// option parameters
	static int numberOfExercisePeriods = 20;
	static double optionPeriodLength = 1;
	private static double swaprate = 0.02;

	@Test
	public void testWithManyPaths() throws CalculationException {
		main(null);
	}

	public static void main(String[] args) throws CalculationException {
		// set parameters
		CreateTestModel.setLastTime(lastTimePoint);
		CreateTestModel.setLiborRateTimeHorzion(lastTimePoint);
		CreateTestModel.setTimeDiscretizationPeriodLength(timeDiscretizationLength);
		CreateTestModel.setLiborPeriodLength(liborPeriodLength);

		CreateTestModel.setNumberOfPaths(numberOfPaths);

		TestValuationMethods.setNumberOfSubsimulationsStepA(numberOfSubsimulationsStepA);
		TestValuationMethods.setNumberOfSubsimulationsStepB(numberOfSubsimulationsStepB);
		TestValuationMethods.setTolerance((double) 1 / (double) numberOfPaths);

		TestValuationMethods.setNumberOfExercisePeriods(numberOfExercisePeriods);
		TestValuationMethods.setOptionPeriodLength(optionPeriodLength);
		TestValuationMethods.setSwaprate(swaprate );

		System.out.println("Number of paths: " + numberOfPaths);
		System.out.println("Number of subsimulation paths step A: " + numberOfSubsimulationsStepA);
		System.out.println("Number of subsimulation paths step B: " + numberOfSubsimulationsStepB);
		System.out.println("Number of exercise periods: " + numberOfExercisePeriods);

		System.out.println("Time discretization period length: " + timeDiscretizationLength);
		System.out.println("LIBOR period length: " + liborPeriodLength);
		System.out.println("Option period length: " + optionPeriodLength);
		System.out.println("Option period length: " + swaprate);
		

		// run test

		TestValuationMethods.executePrintSwaptionValuationMethods();

	}

}
