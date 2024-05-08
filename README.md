# Weight loss analysis project

This repository contains code for analyzing weight loss data in mice models in response to multiple type of infections. The analysis includes statistical tests, mixed-effect modeling, simulations, and power calculations. The project aims to investigate the optimal weight loss threshold to prevent excessive suffering in mice while maintaining statistical power.


For more details, please refer to the article [here](link_to_article).

## Files description

1. `statistical_analysis.ipynb`: Jupyter Notebook containing code for standard statistical analysis on the data.
2. `statistical_analysis_over_time.ipynb`: Jupyter Notebook for analyzing weight loss over time.
3. `mixed_effect_model.ipynb`: Jupyter Notebook implementing mixed-effect modeling (comparing survivor mice to non survivor mice).
4. `run_simulation.ipynb`: Jupyter Notebook for running simulations in order to compute the statistcal power applied to different thresold of sacrifice.
5. `analyse_power_in_house.py`: Python script for conducting in-house power calculations, i.e re-evaluate the statistical test (log rank) in each experiments with new threshold.
6. `determine_mouse_death.py`: Python script to compute at which time a mouse is dead with new threshold of sacrifice.
7. `simulation.py`: Python script implementing the simulations functions.
8. `graphiques.ipynb`: Jupyter Notebook containing code for generating graphs and visualizations of for the article.

## Usage

To run the analysis:

1. Clone the repository to your local machine.
2. Ensure you have Python and Jupyter Notebook installed.
3. Open the desired Jupyter Notebook files in your preferred environment.
4. Follow the instructions within each notebook to execute the analysis steps.

## Contributors

- Maëlick Brochut

## License

This project is licensed under the [License Name] License - see the [LICENSE](LICENSE) file for details.
