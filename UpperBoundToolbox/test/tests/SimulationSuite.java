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
	static int numberOfPaths = 10000;
	static int numberOfSubsimulationsStepA = 10000;
	static int numberOfSubsimulationsStepB = 10000;
	// option parameters
	static int numberOfExercisePeriods = 20;
	static double optionPeriodLength=1;

	@Test
	public void testWithManyPaths() throws CalculationException {
		main(null);
	}

	public static void main(String[] args) throws CalculationException {

		double swapRate = 0.02;
		// out of the money option
		double rateShift = -0.005;
		testOptionFlatCurve(rateShift,swapRate);
		// at the money option
		rateShift=0;
		testOptionFlatCurve(rateShift,swapRate);
		// in  the money option
		rateShift = +0.005;
		testOptionFlatCurve(rateShift,swapRate);	
	}

	private static void testOptionFlatCurve(double rateShift, double swapRate) throws CalculationException
	{
	double[] forwardInterpolationRates = new double[]{swapRate-rateShift,swapRate-rateShift};
	double[] forwardInterpolationTimePoints = new double[] {0,lastTimePoint};

	executeSimulationRunWithFlexibleRates(swapRate,forwardInterpolationRates,forwardInterpolationTimePoints);
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
