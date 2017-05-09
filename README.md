# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

---
![Lidar and Radar Filter Visual](/media/Position_Velocity_RL.png)

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
   - On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../out_data/out.txt`

## Rubric Points
### Compiling
The code can be compiled using the basic build instructions above.  Here is an example of the build on my machine:
![build](/media/build.png)



### Accuracy
The table below shows the RMSE values for the filter when using the combinations of Lidar and Radar data.  The rubric RMSE values are shown as a maximum acceptable.


| Rubric RMSE | Lidar & Radar | Lidar  | Radar  |
| :---------: | :-----------: | :----: | :----: |
|   0.09 >=   |    0.0699     | 0.1232 | 0.1541 |
|   0.10 >=   |    0.0810     | 0.0977 | 0.2137 |
|   0.40 >=   |    0.2878     | 0.6733 | 0.3218 |
|   0.30 >=   |    0.1873     | 0.2420 | 0.2560 |



### Algorithm

#### First Measurement

The state and covariance are initialized based on the first measurement.  For Radar the measurement is converted into the state space and initialized.  The Lidar measurement is initialized directly.  Lidar and Radar have unique initial covariance matrices due to their unique measurement variance.  Initialization is done in `ukf.cpp` in the `ProcessMeasurement` on lines 88-132.

#### Prediction / Update

The prediction / update cycle is defined in `ukf.cpp` in the `ProcessMeasurement` function on line 88.  After the initialization is completed on the first measurement, `ProcessMeasurement1` first runs `Prediction` on line 142 then runs either `UpdateRadar` or `UpdateLidar`  update based on the measurement type on lines 147 and 152 respectively.

#### Radar / Lidar

It is possible to handle Radar or Lidar data in the same data stream.  `UpdateRadar` requires the extra step of conversion to the measurement space while `UpdateLidar` does not.  There are enough differences between the two update procedures that it was decided to keep them separate to reduce the number of control flow statements.

### Code Efficiency

Efforts were made to minimize redundant calculations and control flow statements.  Currently the filter is able to process 379 measurements per second.

## Visualize Results

The filter allows for processing any combination of the Lidar and Radar measurements in the data stream.  Setting what sensor measurements to process is done using the flags in `ukf.cpp`:

16: `use_laser_ = true;`

19: ` use_radar_ = true;`

Here are the visual results of the combinations:

### Lidar and Radar

![Lidar and Radar Filter Visual](/media/Position_Velocity_RL.png)

### Lidar Only

![Lidar and Radar Filter Visual](/media/Position_Velocity_L.png)

### Radar Only

![Lidar and Radar Filter Visual](/media/Position_Velocity_R.png)



## NIS

NIS or Normalized Innovation Squared is a quick way to understand if we are over or under estimating the uncertainty in our system.  NIS was used to tune the linear acceleration, `std_a_` ,and yaw acceleration, `std_yawdd_`, terms.  The Lidar and Radar measurements were processed separately and the process noise parameters were adjusted such that about only 5% of NIS values were observed over the 95% X2 line.

### Radar NIS

![Radar NIS](/media/NIS-R.png)

Here we can see that the Lidar NIS is 0 since we did not process Lidar measurements in this analysis.  The Radar measurement is showing good conformance to the expected X2 distribution with 5.2% above the 5% X2 line and 6.0% below the 95% X2 line.

### Laser NIS

![Radar NIS](/media/NIS-L.png)

Here we can see that the Radar NIS is 0 since we did not process Radar measurements in this analysis.  The Laser measurement is showing good conformance to the expected X2 distribution with 6.0% above the 5% X2 line and 7.2% below the 95% X2 line.  Further adjustment of the noise parameters may improve this but the filter is producing excellent results as is.

