package bermudanSwaptionFramework;

import net.finmath.stochastic.RandomVariable;

/**
 * Interface to enable the use of several exercise strategies.
 * @author Lennart Quante
 * @version 1.0
 */
public interface ExerciseStrategyInterface {

	/**
	 * Return the exercise indicator of a given Bermudan option using a given
	 * continuation and exercise value. The definition of the exercise strategy
	 * method is part of the object implementing this interface.
	 *
	 * @param continuationValue Given continuation value of the Bermudan option.
	 * @param exerciseValue     Given exercise value of the Bermudan option.
	 * 
	 * @return The random variable of exercise indicators
	 */

	RandomVariable getExerciseIndicators(RandomVariable continuationValue, RandomVariable exerciseValue);
}
