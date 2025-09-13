// Sapien Robotics Q1: Multi-Wheel Odometry with Fixed-Point Kalman (ESP32 Arduino)
// Author: Arun Kumar Gupta, Date: September 13, 2025

#include <Arduino.h>

// Pin Definitions (Hall Sensors)
#define FL_PIN 25  // Front Left
#define RL_PIN 26  // Rear Left
#define FR_PIN 27  // Front Right
#define RR_PIN 33  // Rear Right

// Parameters
#define DEBOUNCE_MS 5
#define WHEEL_RADIUS_MM 50    // mm
#define BASE_WIDTH_MM 300     // mm
#define PULSES_PER_REV 20
#define DT_MS 20              // 50Hz
#define FIXED_Q 15
#define FIXED_SCALE (1L << FIXED_Q)  // 32768
#define TO_FIXED(x) ((int32_t)((x) * FIXED_SCALE))
#define FROM_FIXED(x) ((float)(x) / FIXED_SCALE)
#define PI_FIXED TO_FIXED(3.14159265359)

// Volatile Counters and Timestamps
volatile int32_t fl_count = 0, rl_count = 0, fr_count = 0, rr_count = 0;
volatile unsigned long last_fl = 0, last_rl = 0, last_fr = 0, last_rr = 0;

// State Vector [x, y, theta, v, omega]
int32_t state[5] = {0};  // mm, mm, rad, mm/s, rad/s

// Covariance Matrix P (5x5), Process Noise Q, Measurement Noise R
int32_t P[5][5] = {
  {TO_FIXED(10), 0, 0, 0, 0},
  {0, TO_FIXED(10), 0, 0, 0},
  {0, 0, TO_FIXED(0.1), 0, 0},
  {0, 0, 0, TO_FIXED(100), 0},
  {0, 0, 0, 0, TO_FIXED(0.1)}
};
int32_t Q[5] = {TO_FIXED(0.01), TO_FIXED(0.01), TO_FIXED(0.05), TO_FIXED(0.1), TO_FIXED(0.1)};
int32_t R[2] = {TO_FIXED(0.05), TO_FIXED(0.05)};

unsigned long last_out = 0;
int32_t dt_fixed = TO_FIXED(DT_MS / 1000.0);

// Sin/Cos Lookup Table (Q15, 0 to PI/2, 64 entries)
const int32_t sin_lut[65] = {
  0, 804, 1608, 2411, 3213, 4013, 4812, 5609, 6405, 7198, 7989, 8778, 9564, 10347, 11127, 11805,
  12480, 13152, 13821, 14487, 15150, 15810, 16467, 17121, 17772, 18420, 19065, 19707, 20346, 20982,
  21615, 22244, 22871, 23494, 24114, 24731, 25345, 25956, 26563, 27167, 27768, 28365, 28959, 29549,
  30136, 30720, 31300, 31877, 32450, 33020, 33586, 34149, 34708, 35264, 35815, 36363, 36907, 37448,
  37985, 38519, 39049, 39576, 40099, 40619, 41135, 41648, 42156, 42662, 43164, 43663, 44158, 44650,
  45138, 45623, 46104, 46582, 47056, 47527, 47994, 48458, 48918, 49375, 49828, 50278, 50724, 51167,
  51606, 52042, 52474, 52903, 53328, 53750, 54168, 54583, 54994, 55402, 55806, 56207, 56604, 56998,
  57388, 57775, 58158, 58538, 58914, 59287, 59656, 60022, 60384, 60743, 61098, 61450, 61798, 62143,
  62484, 62822, 63156, 63487, 63814, 64138, 64458, 64775, 65088, 65398, 65604, 65807, 66006, 66199,
  66390, 66576, 66759, 66938, 67114, 67286, 67455, 67620, 67782, 67940, 68094, 68245, 68392, 68536,
  68676, 68813, 68946, 69076, 69102, 69124, FIXED_SCALE
};

// Helper Functions
int32_t fp_mul(int32_t a, int32_t b) { return (int64_t)a * b >> FIXED_Q; }
int32_t fp_add(int32_t a, int32_t b) { return a + b; }
int32_t fp_sin(int32_t rad) {
  int32_t norm = (rad % (2 * PI_FIXED) + 2 * PI_FIXED) % (2 * PI_FIXED);
  if (norm > PI_FIXED) norm = 2 * PI_FIXED - norm;
  if (norm > (PI_FIXED / 2)) norm = PI_FIXED - norm;
  int32_t idx = fp_mul(norm, TO_FIXED(128.0 / 3.14159265359));
  int idx_int = idx >> FIXED_Q;
  if (idx_int > 63) idx_int = 63;
  int32_t frac = idx & ((1 << FIXED_Q) - 1);
  int32_t lo = sin_lut[idx_int];
  int32_t hi = sin_lut[idx_int + 1];
  return lo + fp_mul((hi - lo), frac >> (FIXED_Q - 6));
}
int32_t fp_cos(int32_t rad) { return fp_sin(rad + (PI_FIXED / 2)); }

// 2x2 Matrix Inverse
void mat_inv2(int32_t m[2][2], int32_t inv[2][2]) {
  int64_t det = (int64_t)m[0][0] * m[1][1] - (int64_t)m[0][1] * m[1][0];
  if (det == 0) return;
  int32_t idet = FIXED_SCALE * FIXED_SCALE / (det >> FIXED_Q);
  inv[0][0] = fp_mul(m[1][1], idet);
  inv[1][1] = fp_mul(m[0][0], idet);
  inv[0][1] = -fp_mul(m[0][1], idet);
  inv[1][0] = -fp_mul(m[1][0], idet);
}

// Interrupt Service Routines (ISRs)
void IRAM_ATTR fl_isr() {
  unsigned long now = millis();
  if (now - last_fl > DEBOUNCE_MS) { fl_count++; last_fl = now; }
}
void IRAM_ATTR rl_isr() {
  unsigned long now = millis();
  if (now - last_rl > DEBOUNCE_MS) { rl_count++; last_rl = now; }
}
void IRAM_ATTR fr_isr() {
  unsigned long now = millis();
  if (now - last_fr > DEBOUNCE_MS) { fr_count++; last_fr = now; }
}
void IRAM_ATTR rr_isr() {
  unsigned long now = millis();
  if (now - last_rr > DEBOUNCE_MS) { rr_count++; last_rr = now; }
}

// Kalman Predict
void kalman_predict() {
  int32_t c = fp_cos(state[2]);
  int32_t s = fp_sin(state[2]);
  int32_t delta_d = fp_mul(state[3], dt_fixed);  // v * dt
  state[0] = fp_add(state[0], fp_mul(delta_d, c));
  state[1] = fp_add(state[1], fp_mul(delta_d, s));
  state[2] = fp_add(state[2], fp_mul(state[4], dt_fixed));  // theta += omega * dt
  // Propagate P (simplified: add Q)
  for (int i = 0; i < 5; i++) P[i][i] += Q[i];
}

// Kalman Update
void kalman_update(int32_t dist_l, int32_t dist_r) {
  int32_t half_base = TO_FIXED(BASE_WIDTH_MM / 2000.0);  // mm to m scale
  int32_t H[2][5] = {{0, 0, 0, dt_fixed, fp_mul(-dt_fixed, half_base)},
                     {0, 0, 0, dt_fixed, fp_mul(dt_fixed, half_base)}};

  // Predicted measurement
  int32_t pred[2] = {fp_mul(state[3], dt_fixed) + fp_mul(state[4], fp_mul(-dt_fixed, half_base)),
                     fp_mul(state[3], dt_fixed) + fp_mul(state[4], fp_mul(dt_fixed, half_base))};

  // Innovation
  int32_t y[2] = {dist_l - pred[0], dist_r - pred[1]};

  // Innovation covariance S (approx 2x2)
  int32_t S[2][2] = {{P[3][3], P[3][4]}, {P[4][3], P[4][4]}};
  S[0][0] += R[0]; S[1][1] += R[1];

  // Inverse of S
  int32_t S_inv[2][2];
  mat_inv2(S, S_inv);

  // Kalman Gain K (5x2, approx for v, omega)
  int32_t K[5][2] = {{0}};
  K[3][0] = fp_mul(P[3][3], S_inv[0][0]) + fp_mul(P[3][4], S_inv[1][0]);
  K[3][1] = fp_mul(P[3][3], S_inv[0][1]) + fp_mul(P[3][4], S_inv[1][1]);
  K[4][0] = fp_mul(P[4][3], S_inv[0][0]) + fp_mul(P[4][4], S_inv[1][0]);
  K[4][1] = fp_mul(P[4][3], S_inv[0][1]) + fp_mul(P[4][4], S_inv[1][1]);

  // Update state (v, omega mainly)
  state[3] += fp_mul(K[3][0], y[0]) + fp_mul(K[3][1], y[1]);
  state[4] += fp_mul(K[4][0], y[0]) + fp_mul(K[4][1], y[1]);

  // Update P (simplified Joseph form)
  for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) P[i][j] = (P[i][j] * 900LL) / 1000;
  for (int i = 0; i < 5; i++) P[i][i] += Q[i];
}

void setup() {
  Serial.begin(115200);
  pinMode(FL_PIN, INPUT_PULLUP); attachInterrupt(digitalPinToInterrupt(FL_PIN), fl_isr, RISING);
  pinMode(RL_PIN, INPUT_PULLUP); attachInterrupt(digitalPinToInterrupt(RL_PIN), rl_isr, RISING);
  pinMode(FR_PIN, INPUT_PULLUP); attachInterrupt(digitalPinToInterrupt(FR_PIN), fr_isr, RISING);
  pinMode(RR_PIN, INPUT_PULLUP); attachInterrupt(digitalPinToInterrupt(RR_PIN), rr_isr, RISING);
  Serial.println("Odometry System Ready");
}

void loop() {
  unsigned long now = millis();
  if (now - last_out >= DT_MS) {
    // Atomic read and reset counts
    noInterrupts();
    int32_t dfl = fl_count; fl_count = 0;
    int32_t drl = rl_count; rl_count = 0;
    int32_t dfr = fr_count; fr_count = 0;
    int32_t drr = rr_count; rr_count = 0;
    interrupts();

    // Fuse and calculate distances (mm)
    int32_t delta_l = (dfl + drl) / 2;
    int32_t delta_r = (dfr + drr) / 2;
    int32_t dist_l = fp_mul(delta_l, TO_FIXED(2 * 3.14159 * WHEEL_RADIUS_MM / PULSES_PER_REV));
    int32_t dist_r = fp_mul(delta_r, TO_FIXED(2 * 3.14159 * WHEEL_RADIUS_MM / PULSES_PER_REV));

    // Kalman Filter
    kalman_predict();
    kalman_update(dist_l, dist_r);

    // Output (convert to m, rad, m/s, rad/s)
    Serial.printf("%.3f:%.3f:%.3f:%.3f:%.3f\n",
                  FROM_FIXED(state[0]) / 1000.0, FROM_FIXED(state[1]) / 1000.0,
                  FROM_FIXED(state[2]), FROM_FIXED(state[3]) / 1000.0, FROM_FIXED(state[4]));

    last_out = now;
  }
}
