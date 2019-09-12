package simulationMethods;

import net.finmath.exception.CalculationException;

/**
 * Main class to execute simulations
 * @author Lennart Quante
 *
 */
public class ExecuteSimulations {

	// model parameters
	static double lastTimePoint = 30;
	static double timeDiscretizationLength = 0.25;
	static double liborPeriodLength = 0.5;
	// monte carlo parameters
	static int numberOfPaths = 10000;
	static int numberOfSubsimulationsStepA = 1000;
	static int numberOfSubsimulationsStepB = 1000;
	// option parameters
	static int numberOfExercisePeriods = 20;
	static double optionPeriodLength = 1;

	public static void main(String[] args) throws CalculationException {
		SimulationFactory testSimulation = new SimulationFactory(lastTimePoint, timeDiscretizationLength, liborPeriodLength,
				numberOfPaths, numberOfSubsimulationsStepA, numberOfSubsimulationsStepB, numberOfExercisePeriods,
				optionPeriodLength);
		double swapRate = 0.02;
		// out of the money option
		double rateShift = -0.005;
		testSimulation.testOptionFlatCurve(rateShift, swapRate);
		// at the money option
		rateShift = 0;
		testSimulation.testOptionFlatCurve(rateShift, swapRate);
		// in the money option
		rateShift = +0.005;
		testSimulation.testOptionFlatCurve(rateShift, swapRate);
	}

}
