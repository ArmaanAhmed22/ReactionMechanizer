"""Contains tools to visualize a reaction.
"""
from enum import Enum
from typing import Any, Dict, List, Tuple, Union
import typing
from reaction_mechanizer.pathway.reaction import DifferentialEquationModel, ReactionMechanism, SimpleStep
from scipy.integrate import odeint, solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
import pandas as pd


class ReactionEvent(Enum):
    """Enum for the possible reaction events, such as:

    CHANGE_CONCENTRATION: how much the concentration of the species should be changed by.
        Additional Info: Tuple(str: species of interest, float: change in concentration)

    SET_CONCENTRATION: what to set the concentration of the species to.
        Additional Info: Tuple(str: species of interest, float: new concentration)

    SMOOTH_CHANGE_CONCENTRATION: TBD
    """
    CHANGE_CONCENTRATION = 0
    SET_CONCENTRATION = 1
    SMOOTH_CHANGE_CONCENTRATION = 2


class ReactionVisualizer:
    """Visualier for either `SimpleStep` or `ReactionMechanism`
    """
    def __init__(self, reaction: Union[SimpleStep, ReactionMechanism]):
        """Create a `ReactionVisualizer` object to model the progression of a reaction

        Args:
            reaction (Union[SimpleStep, ReactionMechanism]): The reaction object to model
        """
        self.reaction: Union[SimpleStep, ReactionMechanism] = reaction

    def get_states(self,
                   initial_state: Dict[str, float],
                   time_end: float,
                   number_steps: int,
                   initial_time: float = 0,
                   ode_override: Union[Dict[str, DifferentialEquationModel], None] = None) -> Any:
        """Get concentration of the species in this reaction, with model specifications given.

        Args:
            initial_state (Dict[str, float]): initial concentration of all species in reaction
            time_end (float): The end time for this model
            number_steps (int): The granularity of this model. The higher the number of steps, the more accurate the model.
            initial_time (float, optional): The time to start the model at. Defaults to 0.
            ode_override (Union[Dict[str, DifferentialEquationModel], None], optional):
                Dictionary containing the species to override the differential equation of using the provided one. Defaults to None.

        Returns:
            Any: 2D array where the rows represent the concentrations of the species at different times
                (between `initial_time` and `end_time` and using `number_steps`). The columns are the species in the order given by `initial_state`
        """
        ode_dict: Dict[str, DifferentialEquationModel] = self.reaction.get_differential_equations()
        ode_dict_temp = {key: ode_dict[key] for key in initial_state.keys()}
        ode_dict = ode_dict_temp
        if ode_override is not None:
            for key, ode in ode_override.items():
                ode_dict[key] = ode
        ode_function = _get_simple_step_ode_function(ode_dict, list(initial_state.keys()))

        times = np.linspace(initial_time, time_end, number_steps)

        cur_state = list(initial_state.values())

        return odeint(ode_function, cur_state, times)

    def progress_reaction(self,
                          initial_state: Dict[str, float],
                          time_end: float,
                          number_steps: int,
                          events: Union[List[Tuple[float, ReactionEvent, Tuple[Any]]], None] = None,
                          out: Union[str, None] = None,
                          show_intermediates: bool = True) -> pd.DataFrame:
        """Generate model for reaction

        Args:
            initial_state (Dict[str, float]): initial concentration of all species in reaction
            time_end (float): The end time for this model
            number_steps (int): The granularity of this model. The higher the number of steps, the more accurate the model.
            events (List[Tuple[float, ReactionEvent, Tuple[Any]]], optional):
                The list of events to occur during a specified time in the reaction. \
                    A single event is represented by a tuple holding the time of the perturbation, the type of perturbation (`ReactionEvent`), \
                        and the additional information associated with the `ReactionEvent` selected. Defaults to None.
            out (Union[str, None], optional):
                If a string is added, a png visually representing the reaction is created at the specified location (and a `DataFrame` is returned). \
                    Otherwise, just the `DataFrame` is returned. Defaults to None.
            show_intermediates (bool, optional): Whether to show the intermediate species in the graph (doesn't affect dataframe or `SimpleStep`s).

        Returns:
            pd.DataFrame: DataFrame representing the concentrations of the species in the reaction
        """
        data: Any = np.ndarray((0, 0))
        if events is not None:
            sorted_events = (*sorted(events, key=lambda x: x[0]), (time_end, None, tuple()))
            cur_state = dict(initial_state)
            prev_time_point_discretized: float = 0

            for time_point, reaction_event_type, additional_info in sorted_events:
                # First, discretize the time points so that everything is measured to the granularity of "number_steps"
                time_interval = time_end / number_steps
                time_point_discretized = round(time_point / time_interval) * time_interval

                cur_number_steps = round((time_point - prev_time_point_discretized) / time_end * number_steps)

                cur_data = self.get_states(cur_state, time_point, cur_number_steps, initial_time=prev_time_point_discretized)
                data = cur_data if len(data) == 0 else typing.cast(Any, np.concatenate([data, cur_data]))

                if reaction_event_type == ReactionEvent.CHANGE_CONCENTRATION:
                    for index, key in enumerate(cur_state.keys()):
                        if key == additional_info[0]:
                            cur_state[key] = cur_data[-1, index] + additional_info[1]
                        else:
                            cur_state[key] = cur_data[-1, index]
                elif reaction_event_type == ReactionEvent.SET_CONCENTRATION:
                    for index, key in enumerate(cur_state.keys()):
                        if key == additional_info[0]:
                            cur_state[key] = cur_data[-1, index] + additional_info[1]
                        else:
                            cur_state[key] = cur_data[-1, index]
                elif reaction_event_type == ReactionEvent.SMOOTH_CHANGE_CONCENTRATION:
                    pass
                prev_time_point_discretized = time_point_discretized
        else:
            data = self.get_states(initial_state, time_end, number_steps)
        times = np.linspace(0, time_end, number_steps)

        _, ax = plt.subplots()
        plt.tight_layout()
        if out:
            dont_show: List[str] = []
            if not show_intermediates and type(self.reaction) == ReactionMechanism:
                dont_show.extend(self.reaction.get_intermediates())
            for i, thing in enumerate(initial_state.keys()):
                if thing not in dont_show:
                    sns.lineplot(x=times, y=data[:, i], label="$"+thing+"$", ax=ax)
            ax.legend()
            sns.despine(ax=ax)
            ax.margins(x=0, y=0)
            _, top = ax.get_ylim()
            ax.set_ylim([0, top*1.05])
            plt.savefig(str(out), bbox_inches="tight", dpi=600)
        return pd.DataFrame({"Time": times, **{thing: data[:, i] for i, thing in enumerate(initial_state.keys())}})

    def get_states_robust(self,
                          initial_state: Dict[str, float],
                          time_end: float,
                          initial_time: float = 0,
                          ode_override: Union[Dict[str, DifferentialEquationModel], None] = None,
                          **kwargs_solve_ivp) -> Any:
        """Get concentration of the species in this reaction, with model specifications given.
        "get_states_robust" uses the more diverse solve_ivp numerical solver rather than odeint in "get_states"

        Args:
            initial_state (Dict[str, float]): initial concentration of all species in reaction
            time_end (float): The end time for this model
            number_steps (int): The granularity of this model. The higher the number of steps, the more accurate the model.
            initial_time (float, optional): The time to start the model at. Defaults to 0.
            ode_override (Union[Dict[str, DifferentialEquationModel], None], optional):
                Dictionary containing the species to override the differential equation of using the provided one. Defaults to None.

        Returns:
            Any: 2D array where the rows represent the concentrations of the species at different times
                (between `initial_time` and `end_time` and using `number_steps`). The columns are the species in the order given by `initial_state`
        """
        ode_dict: Dict[str, DifferentialEquationModel] = self.reaction.get_differential_equations()
        ode_dict_temp = {key: ode_dict[key] for key in initial_state.keys()}
        ode_dict = ode_dict_temp
        if ode_override is not None:
            for key, ode in ode_override.items():
                ode_dict[key] = ode
        ode_function = _get_simple_step_ode_function(ode_dict, list(initial_state.keys()), inverse_cur_state_and_t=True)

        cur_state = list(initial_state.values())
        return solve_ivp(ode_function, (initial_time, time_end), cur_state, **kwargs_solve_ivp)

    def progress_reaction_robust(self,
                                 initial_state: Dict[str, float],
                                 time_end: float,
                                 events: Union[List[Tuple[float, ReactionEvent, Tuple[Any]]], None] = None,
                                 out: Union[str, None] = None,
                                 show_intermediates: bool = True,
                                 **kwargs_solve_ivp) -> pd.DataFrame:
        """Generate model for reaction

        Args:
            initial_state (Dict[str, float]): initial concentration of all species in reaction
            time_end (float): The end time for this model
            number_steps (int): The granularity of this model. The higher the number of steps, the more accurate the model.
            events (Union[List[Tuple[float, ReactionEvent, Tuple[Any]]], None], optional):
                The list of events to occur during a specified time in the reaction. \
                    A single event is represented by a tuple holding the time of the perturbation, the type of perturbation (`ReactionEvent`), \
                        and the additional information associated with the `ReactionEvent` selected. Defaults to None.
            out (Union[str, None], optional):
                If a string is added, a png visually representing the reaction is created at the specified location (and a `DataFrame` is returned). \
                    Otherwise, just the `DataFrame` is returned. Defaults to None.
            show_intermediates (bool, optional): Whether to show the intermediate species in the graph (doesn't affect dataframe or `SimpleStep`s).

        Returns:
            pd.DataFrame: DataFrame representing the concentrations of the species in the reaction
        """
        data: Any = np.ndarray((0, 0))
        if events is not None:
            print(events)
            sorted_events = (*sorted(events, key=lambda x: x[0]), (time_end, None, tuple()))
            cur_state = dict(initial_state)
            prev_time_point: float = 0

            for time_point, reaction_event_type, additional_info in sorted_events:

                cur_data = self.get_states_robust(cur_state, time_point, initial_time=prev_time_point, **kwargs_solve_ivp)
                data = cur_data if len(data) == 0 else typing.cast(Any, np.concatenate([data, cur_data]))

                if reaction_event_type == ReactionEvent.CHANGE_CONCENTRATION:
                    for index, key in enumerate(cur_state.keys()):
                        if key == additional_info[0]:
                            cur_state[key] = cur_data[-1, index] + additional_info[1]
                        else:
                            cur_state[key] = cur_data[-1, index]
                elif reaction_event_type == ReactionEvent.SET_CONCENTRATION:
                    for index, key in enumerate(cur_state.keys()):
                        if key == additional_info[0]:
                            cur_state[key] = cur_data[-1, index] + additional_info[1]
                        else:
                            cur_state[key] = cur_data[-1, index]
                elif reaction_event_type == ReactionEvent.SMOOTH_CHANGE_CONCENTRATION:
                    pass
                prev_time_point = time_point
        else:
            data = self.get_states_robust(initial_state, time_end, **kwargs_solve_ivp)

        _, ax = plt.subplots()
        plt.tight_layout()
        if out:
            dont_show: List[str] = []
            if not show_intermediates and type(self.reaction) == ReactionMechanism:
                dont_show.extend(self.reaction.get_intermediates())
            for i, thing in enumerate(initial_state.keys()):
                if thing not in dont_show:
                    sns.lineplot(x=data.t, y=data.y[i, :], label="$"+thing+"$", ax=ax)
            ax.legend()
            sns.despine(ax=ax)
            ax.margins(x=0, y=0)
            _, top = ax.get_ylim()
            ax.set_ylim([0, top*1.05])
            plt.savefig(str(out), bbox_inches="tight", dpi=600)
        return pd.DataFrame({"Time": data.t, **{thing: data.y[i, :] for i, thing in enumerate(initial_state.keys())}})

    def animate_progress_reaction(self,
                                  video_destination_no_extension: str,
                                  video_length: float,
                                  fps: int,
                                  extension: str = "mp4",
                                  uses_robust: bool = False,
                                  **progress_reaction_args):
        """Generate video visualization for reaction

        Args:
            video_destination_no_extension (str): The path to store the video at
            video_length (float): The length of the resulting video
            fps (int): The fps of the video
            extension (str, optional): The extensions of the video [example values: "mp4", "gif", "mov", etc]. Defaults to "mp4".
            uses_robust (bool, optional): Whether or not to use "self.progress_reaction" or "self.progress_reaction_robust"
        """
        df = self.progress_reaction(**progress_reaction_args) if not uses_robust else self.progress_reaction_robust(**progress_reaction_args)
        if not progress_reaction_args.get("show_intermediates", True) and type(self.reaction) == ReactionMechanism:
            df = df.drop(self.reaction.get_intermediates(), axis="columns")
        writer = animation.writers["ffmpeg"](fps=fps, metadata={"artist": "ReactionMechanizer"}, bitrate=1800)  # Non-python dependency!
        _, dummy_ax = plt.subplots()
        plt.tight_layout()
        df2 = df.drop("Time", axis="columns").stack().reset_index()  # level_1: species, 0: concentration values
        df2 = df2.rename({"level_1": "Species", 0: "Concentration"}, axis="columns")
        df2["Time"] = [cur_time for cur_time in df["Time"] for _ in df2["Species"].unique()]
        sns.lineplot(x="Time", y="Concentration", data=df2, hue="Species", ax=dummy_ax)
        _, top_y = dummy_ax.get_ylim()
        _, top_x = dummy_ax.get_xlim()

        new_fig, new_ax = plt.subplots()

        data_1_item = df2.iloc[:len(df2["Species"].unique())]
        cur_data: Dict[str, Dict[str, List[float]]] = \
            {spec: {"x": [data_1_item["Time"].iloc[0]], "y": [data_1_item["Concentration"].iloc[i]]}
             for i, spec in enumerate(df2["Species"].unique())}
        sns.lineplot(x="Time", y="Concentration", data=data_1_item, hue="Species", ax=new_ax)
        new_ax.set_ylim([0, top_y*1.05])
        new_ax.set_xlim([0, top_x*1.05])
        new_ax.margins(x=0, y=0)
        sns.despine(ax=new_ax)

        def animate(frame_index):
            for i, spec in enumerate(df2["Species"].unique()):
                cur_data[spec]["x"].append(df["Time"].iloc[int(frame_index/(video_length*fps)*df.shape[0])])
                cur_data[spec]["y"].append(df[spec].iloc[int(frame_index/(video_length*fps)*df.shape[0])])
                new_ax.get_lines()[i].set_data(cur_data[spec]["x"], cur_data[spec]["y"])
        ani = animation.FuncAnimation(new_fig, animate, frames=video_length*fps, repeat=True)
        ani.save(f"{video_destination_no_extension}.{extension}", writer=writer)


def _get_simple_step_ode_function(differential_equations: Dict[str, DifferentialEquationModel], state_order: List[str], inverse_cur_state_and_t: bool = False):
    active_odes = {key: val.get_lambda() for key, val in differential_equations.items()}

    def simple_step_ode_function(cur_state, t):
        cur_state_order = state_order
        dict_state = {thing: val for thing, val in zip(cur_state_order, cur_state)}
        out = []
        for _, ode in active_odes.items():
            out.append(ode(**dict_state))
        return out

    def simple_step_ode_function_params_inversed(t, cur_state):
        cur_state_order = state_order
        dict_state = {thing: val for thing, val in zip(cur_state_order, cur_state)}
        out = []
        for _, ode in active_odes.items():
            out.append(ode(**dict_state))
        return out
    return simple_step_ode_function if not inverse_cur_state_and_t else simple_step_ode_function_params_inversed
