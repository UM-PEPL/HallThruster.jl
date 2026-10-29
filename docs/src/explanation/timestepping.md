# Timestepping

Each solver step selects a timestep, advances the heavy species and electrons, updates the adaptive fallback counters, and saves a frame if it reaches the next requested save time.

`SimParams.num_save` counts the initial and final frames as well as the frames between them. 
The solver periodically shortens timesteps so that the soler lands exactly on the requested checkpoint times.
These shortened steps may be smaller than `min_dt`.

Internally, timestep values have three distinct roles: `params.simulation.dt` is the requested timestep, `params.cache.dt_limit[]` is the tability limit computed by the heavy-species solver, and `params.dt[]` is the actual timestep used by the physics in the current step.
Adaptive simulations begin with a small startup step. Later steps clamp the
stability limit between `min_dt` and `max_dt` before applying the checkpointconstraint.
After `max_small_steps` consecutive minimum-sized adaptive steps, the solver temporarily uses the configured `dt` for `max_small_steps` steps,then resumes adaptation.
Startup and checkpoint-shortened steps do not count toward this threshold. Checkpoint constraints apply during fallback as well.

Each saved `Frame.dt` records the actual timestep used by the last completed step, rather than the stability limit or the interval between frames. The initial frame has a timestep of zero.

The solver appends a timestamp only when it saves a frame. If a simulation fails early, the returned solution contains only the checkpoints completed before the failure, with matching timestamp and frame counts.
