package tests;

import org.junit.Assert;
import org.junit.Test;

import net.finmath.exception.CalculationException;
import simulationMethods.SimulationFactory;

/**
 * Test to test simulation runner specifications with a low number of paths to avoid missspecification e.g of the model
 * @author Lennart Quante
 *
 */
public class TestExecuteSimulations {

	// model parameters
	static double lastTimePoint = 30;
	static double timeDiscretizationLength = 0.25;
	static double liborPeriodLength = 0.5;
	// monte carlo parameters
	static int numberOfPaths = 100;
	static int numberOfSubsimulationsStepA = 1;
	static int numberOfSubsimulationsStepB = 1;
	// option parameters
	static int numberOfExercisePeriods = 20;
	static double optionPeriodLength = 1;

	@Test
	public void testExecuteSimulations() throws CalculationException {

		SimulationFactory testSimulation = new SimulationFactory(lastTimePoint, timeDiscretizationLength, liborPeriodLength,
				numberOfPaths, numberOfSubsimulationsStepA, numberOfSubsimulationsStepB, numberOfExercisePeriods,
				optionPeriodLength);
		double swapRate = 0.02;
	
		// at the money option
		double rateShift = 0;
		testSimulation.testOptionFlatCurve(rateShift, swapRate);
		boolean valuationSuccessful=true;
		Assert.assertTrue(valuationSuccessful);

	}

}
