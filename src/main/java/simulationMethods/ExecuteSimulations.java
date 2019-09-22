package simulationMethods;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel.Scheme;

/**
 * Main class to execute simulations
 * @author Lennart Quante
 * @version 1.0
 */
public class ExecuteSimulations {

	// model parameters
	static double lastTimePoint = 30;
	static double timeDiscretizationLength = 0.5;
	static double liborPeriodLength = 0.5;
	// option parameters
	static int numberOfExercisePeriods = 20;
	static double optionPeriodLength = 1;
	private static Scheme subsimulationScheme= Scheme.EULER_FUNCTIONAL;

	/**
	 * @param args the monte carlo parameters as string: number of paths, number of Subsimulations A, number of subsimulations B
	 * @throws CalculationException if some rate calculation fails
	 */
	public static void main(String[] args) throws CalculationException {
	
		int numberOfPaths = Integer.parseInt(args[0]);
		int numberOfSubsimulationsStepA = Integer.parseInt(args[1]);
		int numberOfSubsimulationsStepB = Integer.parseInt(args[2]);
		
		SimulationFactory testSimulation = new SimulationFactory(lastTimePoint, timeDiscretizationLength, liborPeriodLength,
				numberOfPaths, numberOfSubsimulationsStepA, numberOfSubsimulationsStepB, numberOfExercisePeriods,
				optionPeriodLength, subsimulationScheme);
		double swapRate = 0.02;
		// out of the money option
		double rateShift = -0.005;
		//testSimulation.testOptionFlatCurve(rateShift, swapRate);
		// at the money option
		rateShift = 0;
		testSimulation.testOptionFlatCurve(rateShift, swapRate);
		// in the money option
		rateShift = +0.005;
		//testSimulation.testOptionFlatCurve(rateShift, swapRate);
	}

}
