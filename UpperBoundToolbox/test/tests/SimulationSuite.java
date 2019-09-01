package tests;

import java.util.Arrays;

import org.junit.Test;

import net.finmath.exception.CalculationException;

public class SimulationSuite {
	// model parameters
	static double lastTimePoint = 30;
	static double timeDiscretizationLength=0.25;
	static double liborPeriodLength = 0.5;
	// monte carlo parameters
	static int numberOfPaths = 10;
	static int numberOfSubsimulationsStepA = 100;
	static int numberOfSubsimulationsStepB = 100;
	// option parameters
	static int numberOfExercisePeriods = 20;
	static double optionPeriodLength=1;

	@Test
	public void testWithManyPaths() throws CalculationException {
		main(null);
	}

	public static void main(String[] args) throws CalculationException {

		// at the money option
		double swaprate = 0.02;
		double [] forwardInterpolationRates= {swaprate,swaprate};
		double [] forwardInterpolationTimePoints= {0,lastTimePoint};

		executeSimulationRunWithFlexibleRates(swaprate,forwardInterpolationRates,forwardInterpolationTimePoints);

		// out of the money option
		swaprate = 0.02;
		forwardInterpolationRates= new double[]{swaprate-0.001,swaprate-0.001};
		forwardInterpolationTimePoints= new double[] {0,lastTimePoint};

		executeSimulationRunWithFlexibleRates(swaprate,forwardInterpolationRates,forwardInterpolationTimePoints);
	}

	public static void executeSimulationRunWithFlexibleRates (double swaprate,double[] forwardInterpolationRates,double[] forwardInterpolationTimePoints) throws CalculationException
	{
		// set parameters
		CreateTestModel.setLastTime(lastTimePoint);
		CreateTestModel.setLiborRateTimeHorzion(lastTimePoint);
		CreateTestModel.setTimeDiscretizationPeriodLength(timeDiscretizationLength);
		CreateTestModel.setLiborPeriodLength(liborPeriodLength);


		CreateTestModel.setForwardInterpolationRates(forwardInterpolationRates);
		CreateTestModel.setForwardInterpolationTimePoints(forwardInterpolationTimePoints);

		CreateTestModel.setNumberOfPaths(numberOfPaths);

		TestValuationMethods.setNumberOfSubsimulationsStepA(numberOfSubsimulationsStepA);
		TestValuationMethods.setNumberOfSubsimulationsStepB(numberOfSubsimulationsStepB);

		CreateTestBermudanSwaption.setNumberOfExercisePeriods(numberOfExercisePeriods);
		CreateTestBermudanSwaption.setPeriodLength(optionPeriodLength);
		CreateTestBermudanSwaption.setSwaprate(swaprate);

		System.out.println("Number of paths: "+numberOfPaths);
		System.out.println("Number of subsimulation paths step A: "+numberOfSubsimulationsStepA);
		System.out.println("Number of subsimulation paths step B: "+numberOfSubsimulationsStepB);
		System.out.println("Number of exercise periods: "+numberOfExercisePeriods);

		System.out.println("Time discretization period length: "+timeDiscretizationLength);
		System.out.println("LIBOR period length: "+liborPeriodLength);
		System.out.println("Initial forward rates: "+Arrays.toString(forwardInterpolationRates));


		System.out.println("Option period length: "+optionPeriodLength);
		System.out.println("Option swap rate: "+swaprate);
		// run test

		TestValuationMethods.executePrintSwaptionValuationMethods();
	}

}
