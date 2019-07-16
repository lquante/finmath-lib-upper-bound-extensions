package drafts;

import net.finmath.stochastic.RandomVariable;

public class SimplestExerciseStrategy implements ExerciseStrategyInterface {

	@Override
	public RandomVariable getExerciseIndicators(RandomVariable continuationValue, RandomVariable exerciseValue) {

		return exerciseValue.sub(continuationValue);
	}

}
