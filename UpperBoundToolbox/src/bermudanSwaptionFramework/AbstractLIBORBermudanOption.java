package bermudanSwaptionFramework;

import java.util.HashMap;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.stochastic.RandomVariable;

/**
 * Provides a basic object for Bermudan interest rate options.
 * @author Lennart Quante
 * @version 1.0
 */
public abstract class AbstractLIBORBermudanOption extends AbstractLIBORMonteCarloProduct {

	private final boolean[] isPeriodStartDateExerciseDate; // Exercise date
	private final double[] fixingDates; // Vector of fixing dates (must be sorted)
	private final double[] periodLengths; // Vector of period lengths (could be ommitted?!)
	private final double[] paymentDates; // Vector of payment dates (same length as fixing dates)
	private final double[] periodNotionals; // Vector of notionals for each period
	private final boolean isCallable; // If true: the right to enter into the underlying product, else the right to
										// terminate the product.

	/**
	 * @param isPeriodStartDateExerciseDate Boolean vector, true if this period is
	 *                                      exercise period
	 * @param fixingDates                   Array of fixing dates
	 * @param periodLengths                 Array of period lengths
	 * @param paymentDates                  Array of payment dates
	 * @param periodNotionals               Array of notionals
	 * @param isCallable                    True if the option is callable,
	 *                                      otherwise cancelable
	 * @param currency The optional currency name of the option as a string.                                     
	 */

	public AbstractLIBORBermudanOption(String currency, boolean[] isPeriodStartDateExerciseDate, double[] fixingDates,
			double[] periodLengths, double[] paymentDates, double[] periodNotionals, boolean isCallable) {
		super(currency);
		this.isPeriodStartDateExerciseDate = isPeriodStartDateExerciseDate;
		this.fixingDates = fixingDates;
		this.periodLengths = periodLengths;
		this.paymentDates = paymentDates;
		this.periodNotionals = periodNotionals;
		this.isCallable = isCallable;
	}

	@Override
	public RandomVariable getValue(double evaluationTime, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {
		
		throw(new CalculationException("The general valuation of Bermudan options is not yet implemented."));
		
		
	}
	

	@Override
	public Map<String, Object> getValues(double evaluationTime, LIBORModelMonteCarloSimulationModel model)
			throws CalculationException {
		RandomVariable value = getValue(evaluationTime, model);
		Map<String, Object> result = new HashMap<>();
		result.put("value", value.getAverage());
		result.put("error", value.getStandardError());
		return result;
	}

	public boolean[] getIsPeriodStartDateExerciseDate() {
		return isPeriodStartDateExerciseDate;
	}

	public double[] getFixingDates() {
		return fixingDates;
	}

	public double[] getPeriodLengths() {
		return periodLengths;
	}

	public double[] getPaymentDates() {
		return paymentDates;
	}

	public double[] getPeriodNotionals() {
		return periodNotionals;
	}

	public boolean isCallable() {
		return isCallable;
	}

}
