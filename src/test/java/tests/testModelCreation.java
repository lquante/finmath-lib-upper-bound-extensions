package tests;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import simulationMethods.TestModelFactory;

public class testModelCreation {

/**	This class tests if a creation of a LMM using the specified factory suceeds:
		 *
		 * @author Lennart Quante
		 */	
	
	@Test
	public void testModel() throws CalculationException {

		TestModelFactory modelFactory = new TestModelFactory();
		LIBORModelMonteCarloSimulationModel liborModel = TestModelFactory.createLIBORMarketModel();

		// print some model properties
		TimeDiscretization timeDiscretization = liborModel.getTimeDiscretization();
		TimeDiscretization liborDiscretization = liborModel.getLiborPeriodDiscretization();
		System.out.println("TimeDiscretization:" + Arrays.toString(timeDiscretization.getAsDoubleArray()));
		System.out.println("LiborDiscretization:" + Arrays.toString(liborDiscretization.getAsDoubleArray()));

		boolean liborCalculated = false;
		// try to calculate all LIBOR rates
		for (int i = 0; i < liborDiscretization.getNumberOfTimeSteps(); i++)
			try {
				RandomVariable[] libors = liborModel.getLIBORs(i);
				System.out.println((i + 1) + "thLIBOR rate " + libors[i].getAverage());
				liborCalculated = true;
			} catch (CalculationException e) {
				liborCalculated = false;
				System.out.println("calculation of libors failed at timepoint" + liborDiscretization.getTime(i));
				e.printStackTrace();

			}

		Assert.assertTrue(liborCalculated);

	}
}
