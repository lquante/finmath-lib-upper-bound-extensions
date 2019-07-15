package drafts;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.stochastic.RandomVariable;

public interface BermudanSwaptionValueEstimatorInterface {

	/**
	 * Return the value estimation of a given Bermudan option using a given
	 * LIBORMonteCarloModel and evaluation time. The definition of the valuation
	 * method is part of the object implementing this interface.
	 *
	 * @param bermudanOption Given Bermudan option.
	 * @param model          The used LIBOR model.
	 * @param evaluationTime The time on which the evaluation should be performed.
	 * @return The value random variable of <code>bermudanOption</code>.
	 * @throws CalculationException
	 */

	RandomVariable getValueEstimation(BermudanSwaption bermudanOption, double evaluationTime,
			LIBORModelMonteCarloSimulationModel model, RandomVariable triggerValues) throws CalculationException;
}
