package bermudanSwaptionFramework;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

/**
 * Interface to enable valuation using different methods through this interface.
 * @author Lennart Quante
 *
 */
public interface BermudanSwaptionValueEstimatorInterface {

	/**
	 * Returns the value estimation of a given Bermudan option using a given
	 * LIBORMonteCarloModel and evaluation time. The definition of the valuation
	 * method is part of the object implementing this interface.
	 *
	 * @param bermudanOption Given Bermudan option.
	 * @param model          The used LIBOR model.
	 * @param evaluationTime The time on which the evaluation should be performed.
	 * @param triggerValues to be used as externaly given
	 * @return The value random variable of <code>bermudanOption</code>.
	 * @throws CalculationException
	 */

	RandomVariable getValueEstimation(BermudanSwaption bermudanOption, double evaluationTime,
			LIBORModelMonteCarloSimulationModel model, RandomVariable triggerValues) throws CalculationException;
}
