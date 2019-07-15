/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 29.02.2008
 */
package modifiedFinmathFilesProcess;

import net.finmath.exception.CalculationException;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;

/**
 * The interface for a stochastic process <i>X</i>.
 *
 * @author Christian Fries
 * @version 1.0
 */
public interface Process {

	/**
	 * This method returns the realization of the process for a given time index.
	 *
	 * @param timeIndex Time index at which the process should be observed
	 * @return The process realizations (given as array of <code>RandomVariable</code>)
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	//	RandomVariable[] getProcessValue(int timeIndex);

	/**
	 * This method returns the realization of a component of the process for a given time index.
	 *
	 * @param timeIndex Time index at which the process should be observed
	 * @param component Component index of the process
	 * @return The process component realizations (given as <code>RandomVariableFromDoubleArray</code>)
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	RandomVariable getProcessValue(int timeIndex, int component) throws CalculationException;

	/**
	 * This method returns the weights of a weighted Monte Carlo method (the probability density).
	 *
	 * @param timeIndex Time index at which the process should be observed
	 * @return A vector of positive weights which sums up to one
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	RandomVariable getMonteCarloWeights(int timeIndex) throws CalculationException;

	/**
	 * @return Returns the numberOfComponents.
	 */
	int getNumberOfComponents();

	/**
	 * @return Returns the timeDiscretizationFromArray.
	 */
	TimeDiscretization getTimeDiscretization();

	/**
	 * @param timeIndex Time index.
	 * @return Returns the time for a given time index.
	 */
	double getTime(int timeIndex);

	/**
	 * Returns the time index for a given simulation time.
	 * @param time The given simulation time.
	 * @return Returns the time index for a given time
	 */
	int getTimeIndex(double time);

	/**
	 * Create and return a clone of this process. The clone is not tied to any model, but has the same
	 * process specification, that is, if the model is the same, it would generate the same paths.
	 *
	 * @return Clone of the process
	 */
	Process clone();

}
