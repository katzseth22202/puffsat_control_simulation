"""Pure C3b terminal-guidance core: ZEM law, tracker noise, plate-frame miss (ADR 0014/0015)."""

from __future__ import annotations

import numpy as np
import pytest

from puffsat_sim.guidance import thrust_command, zem_acceleration


class TestPlateFrameMiss:
    def test_pure_along_track_offset_is_a_time_of_arrival_error_not_a_lateral_miss(
        self,
    ) -> None:
        from puffsat_sim.guidance import plate_frame_miss

        v = (0.0, 10_780.0, 0.0)
        ahead_s = 0.0023
        position = (0.0, 10_780.0 * ahead_s, 0.0)

        miss = plate_frame_miss(
            position, v, toa_s=10_000.0, target_position_m=(0.0, 0.0, 0.0), target_toa_s=10_000.0
        )

        assert miss.lateral_norm_m == pytest.approx(0.0, abs=1e-9)
        assert miss.toa_error_s == pytest.approx(-ahead_s)

    def test_pure_lateral_offset_lands_fully_in_the_plate_plane(self) -> None:
        from puffsat_sim.guidance import plate_frame_miss

        v = (0.0, 10_780.0, 0.0)
        miss = plate_frame_miss(
            (3.0, 0.0, -4.0),
            v,
            toa_s=10_000.2,
            target_position_m=(0.0, 0.0, 0.0),
            target_toa_s=10_000.0,
        )

        assert miss.lateral_norm_m == pytest.approx(5.0)
        assert miss.toa_error_s == pytest.approx(0.2)
        assert np.hypot(*miss.lateral_m) == pytest.approx(5.0)

    def test_mixed_offset_decomposes_exactly(self) -> None:
        from puffsat_sim.guidance import plate_frame_miss

        v = np.array([7_000.0, -5_000.0, 6_000.0])
        unit = v / np.linalg.norm(v)
        lateral = np.cross(unit, [0.0, 0.0, 1.0])
        lateral = 2.5 * lateral / np.linalg.norm(lateral)
        offset = lateral + 40.0 * unit

        miss = plate_frame_miss(
            (offset[0], offset[1], offset[2]),
            (v[0], v[1], v[2]),
            toa_s=0.0,
            target_position_m=(0.0, 0.0, 0.0),
            target_toa_s=0.0,
        )

        assert miss.lateral_norm_m == pytest.approx(2.5)
        assert miss.toa_error_s == pytest.approx(-40.0 / float(np.linalg.norm(v)))


class TestCapture:
    def test_capture_requires_both_the_plate_radius_and_the_toa_window(self) -> None:
        from puffsat_sim.guidance import PlateMiss, capture_fraction

        misses = (
            PlateMiss(lateral_m=(1.0, 1.0), toa_error_s=0.001),
            PlateMiss(lateral_m=(4.0, 4.0), toa_error_s=0.001),
            PlateMiss(lateral_m=(1.0, 0.0), toa_error_s=-0.02),
            PlateMiss(lateral_m=(6.0, 0.0), toa_error_s=0.02),
        )

        assert capture_fraction(misses, plate_radius_m=5.0, toa_limit_s=0.010) == 0.25

    def test_capture_curve_sweeps_the_plate_radius(self) -> None:
        from puffsat_sim.guidance import PlateMiss, capture_curve

        misses = tuple(PlateMiss(lateral_m=(r, 0.0), toa_error_s=0.0) for r in (0.5, 2.0, 4.0, 8.0))

        curve = capture_curve(misses, radii_m=(1.0, 5.0, 10.0))

        assert curve == ((1.0, 0.25), (5.0, 0.75), (10.0, 1.0))


def _run(lateral: tuple[float, float], toa: float = 0.0, dv: float = 0.01) -> object:
    from puffsat_sim.guidance import GuidanceRun, PlateMiss
    from puffsat_sim.terminal import ThrustCommand, executed_plan

    commands = (
        ThrustCommand(start_s=0.0, duration_s=1.0, thrust_n=dv * 25.0, direction=(1.0, 0.0, 0.0)),
    )
    return GuidanceRun(
        miss=PlateMiss(lateral_m=lateral, toa_error_s=toa),
        plan=executed_plan(commands, mass_kg=25.0, saturated=False),
        saturated_fraction=0.0,
    )


class TestGuidanceCell:
    def test_cell_summarizes_rms_and_worst_lateral_and_capture(self) -> None:
        from puffsat_sim.guidance import GuidanceCell

        cell = GuidanceCell(
            label="σ_θ = 10 µrad",
            runs=(_run((3.0, 0.0)), _run((0.0, 4.0)), _run((6.0, 0.0))),
        )

        assert cell.rms_lateral_m == pytest.approx(np.sqrt((9.0 + 16.0 + 36.0) / 3.0))
        assert cell.max_lateral_m == 6.0
        assert cell.capture == pytest.approx(2.0 / 3.0)

    def test_cell_reports_the_worst_run_dv(self) -> None:
        from puffsat_sim.guidance import GuidanceCell

        cell = GuidanceCell(label="x", runs=(_run((0.0, 0.0), dv=0.01), _run((0.0, 0.0), dv=0.04)))

        assert cell.max_dv_m_s == pytest.approx(0.04)


class TestMeasuredCatchRadius:
    def test_largest_entry_offset_still_landing_capture_grade(self) -> None:
        from puffsat_sim.guidance import GuidanceCell, measured_catch_radius_m

        cells = (
            GuidanceCell(label="0 m", runs=(_run((0.01, 0.0)),), axis_value=0.0),
            GuidanceCell(label="300 m", runs=(_run((0.05, 0.0)),), axis_value=300.0),
            GuidanceCell(label="500 m", runs=(_run((1.2, 0.0)),), axis_value=500.0),
            GuidanceCell(label="700 m", runs=(_run((220.0, 0.0)),), axis_value=700.0),
        )

        assert measured_catch_radius_m(cells) == 500.0

    def test_no_cell_within_capture_grade_reads_none(self) -> None:
        from puffsat_sim.guidance import GuidanceCell, measured_catch_radius_m

        cells = (GuidanceCell(label="800 m", runs=(_run((300.0, 0.0)),), axis_value=800.0),)

        assert measured_catch_radius_m(cells) is None


class TestFormatTerminalGuidance:
    def _finding(self) -> object:
        from puffsat_sim.guidance import GuidanceCell, TerminalGuidanceFinding

        return TerminalGuidanceFinding(
            entry_cells=(
                GuidanceCell(label="0 m", runs=(_run((0.002, 0.0)),), axis_value=0.0),
                GuidanceCell(label="400 m", runs=(_run((0.05, 0.0)),), axis_value=400.0),
                GuidanceCell(label="700 m", runs=(_run((250.0, 0.0)),), axis_value=700.0),
            ),
            grade_cells=(
                GuidanceCell(
                    label="σ_θ = 10 µrad",
                    runs=(_run((0.4, 0.0)), _run((0.7, 0.0))),
                    axis_value=10e-6,
                ),
                GuidanceCell(label="σ_rel = 1 m const", runs=(_run((0.3, 0.0)),)),
            ),
            cadence_cells=(GuidanceCell(label="1 Hz", runs=(_run((0.5, 0.0)),), axis_value=1.0),),
            drag_cells=(GuidanceCell(label="Cd ×2.0", runs=(_run((0.6, 0.0)),), axis_value=2.0),),
            cadence_hz=1.0,
            gain=3.0,
            a_max_m_s2=0.016,
            speed_m_s=10_780.0,
        )

    def test_report_reads_the_catch_radius_the_floors_and_the_gates(self) -> None:
        from puffsat_sim.guidance import format_terminal_guidance

        report = format_terminal_guidance(self._finding())

        assert "measured catch radius" in report
        assert "400 m" in report
        assert "σ_θ = 10 µrad" in report
        assert "floor 1.45 m" in report
        assert "capture 100%" in report
        assert "Δv 0.010 m/s" in report
        assert "PASS" in report
        assert "Propellant" in report

    def test_gate_verdicts_pass_at_the_actuator_rails(self) -> None:
        """Riding the rails (thrust at cap, slew recomputed to 1°/s + ε) is not a violation."""
        from puffsat_sim.guidance import (
            GuidanceCell,
            GuidanceRun,
            PlateMiss,
            format_terminal_guidance,
        )
        from puffsat_sim.terminal import FeedforwardPlan, ThrustCommand

        rail_plan = FeedforwardPlan(
            commands=(
                ThrustCommand(start_s=0.0, duration_s=1.0, thrust_n=0.4, direction=(1.0, 0.0, 0.0)),
            ),
            dv_m_s=0.016,
            mass_kg=25.0,
            saturated=True,
            peak_slew_rate_deg_s=1.0 + 1e-12,
        )
        rail_run = GuidanceRun(
            miss=PlateMiss(lateral_m=(0.1, 0.0), toa_error_s=0.0),
            plan=rail_plan,
            saturated_fraction=1.0,
        )
        from dataclasses import replace

        finding = replace(
            self._finding(),
            entry_cells=(GuidanceCell(label="500 m", runs=(rail_run,), axis_value=500.0),),
        )

        report = format_terminal_guidance(finding)

        assert "FAIL" not in report
        assert report.count("PASS") == 2


class TestSlewLimitedDirection:
    def test_within_the_limit_the_commanded_direction_is_taken(self) -> None:
        from puffsat_sim.guidance import slew_limited_direction

        previous = (1.0, 0.0, 0.0)
        commanded = (np.cos(np.radians(0.5)), np.sin(np.radians(0.5)), 0.0)

        result = slew_limited_direction(previous, commanded, max_angle_deg=1.0)

        assert np.allclose(result, commanded)

    def test_beyond_the_limit_the_direction_rotates_by_exactly_the_limit(self) -> None:
        from puffsat_sim.guidance import slew_limited_direction

        previous = (1.0, 0.0, 0.0)
        commanded = (0.0, 1.0, 0.0)

        result = slew_limited_direction(previous, commanded, max_angle_deg=1.0)

        assert np.linalg.norm(result) == pytest.approx(1.0)
        assert np.degrees(np.arccos(np.dot(result, previous))) == pytest.approx(1.0)
        assert result[1] > 0.0

    def test_antiparallel_command_still_rotates_in_a_plane(self) -> None:
        from puffsat_sim.guidance import slew_limited_direction

        result = slew_limited_direction((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), max_angle_deg=2.0)

        assert np.linalg.norm(result) == pytest.approx(1.0)
        assert np.degrees(np.arccos(np.dot(result, (1.0, 0.0, 0.0)))) == pytest.approx(2.0)


class TestHomingFloor:
    def test_matches_the_adr_0015_anchor_at_10_microrad(self) -> None:
        from puffsat_sim.guidance import homing_floor_m

        floor = homing_floor_m(sigma_theta_rad=10e-6, speed_m_s=10_780.0, a_max_m_s2=0.016)

        assert floor == pytest.approx(1.45, abs=0.01)


class TestZemAcceleration:
    def test_commands_k_zem_over_t_go_squared_toward_the_target(self) -> None:
        zem = (300.0, 0.0, -400.0)
        accel = zem_acceleration(zem, t_go_s=100.0, gain=3.0)

        assert np.allclose(accel, (3.0 * 300.0 / 1e4, 0.0, 3.0 * -400.0 / 1e4))


class TestSignificantZem:
    def test_an_estimate_below_the_significance_gate_commands_nothing(self) -> None:
        from puffsat_sim.guidance import significant_zem

        zem = significant_zem((10.0, -20.0, 5.0), knowledge_sigma_m=30.0, t_go_s=200.0)

        assert zem == (0.0, 0.0, 0.0)

    def test_a_significant_estimate_is_acted_on_in_full_not_shrunk(self) -> None:
        from puffsat_sim.guidance import significant_zem

        zem = significant_zem((300.0, 0.0, 400.0), knowledge_sigma_m=30.0, t_go_s=200.0)

        assert zem == (300.0, 0.0, 400.0)

    def test_inside_the_tracking_window_the_gate_is_open_regardless(self) -> None:
        from puffsat_sim.guidance import TRACK_WINDOW_S, significant_zem

        zem = significant_zem((1.0, 0.0, 0.0), knowledge_sigma_m=30.0, t_go_s=TRACK_WINDOW_S - 1.0)

        assert zem == (1.0, 0.0, 0.0)

    def test_perfect_knowledge_is_always_significant(self) -> None:
        from puffsat_sim.guidance import significant_zem

        zem = significant_zem((0.5, 0.0, 0.0), knowledge_sigma_m=0.0, t_go_s=200.0)

        assert zem == (0.5, 0.0, 0.0)


class TestTerminalTick:
    def test_a_significant_zem_fires_the_capped_thrust_along_the_command(self) -> None:
        from puffsat_sim.guidance import terminal_tick

        tick = terminal_tick(
            zem_est_m=(4000.0, 0.0, 0.0),
            knowledge_sigma_m=0.0,
            t_go_s=100.0,
            feedforward_m_s2=(0.0, 0.0, 0.0),
            attitude_dir=None,
            control_period_s=1.0,
            mass_kg=25.0,
        )

        assert tick.fire
        assert tick.saturated
        assert tick.thrust_n == pytest.approx(0.4)
        assert np.allclose(tick.attitude_dir, (1.0, 0.0, 0.0))

    def test_fire_is_held_while_the_gimbal_lags_the_command(self) -> None:
        from puffsat_sim.guidance import terminal_tick

        tick = terminal_tick(
            zem_est_m=(0.0, 4000.0, 0.0),
            knowledge_sigma_m=0.0,
            t_go_s=100.0,
            feedforward_m_s2=(0.0, 0.0, 0.0),
            attitude_dir=(1.0, 0.0, 0.0),
            control_period_s=1.0,
            mass_kg=25.0,
        )

        assert not tick.fire
        angle_deg = np.degrees(np.arccos(np.dot(tick.attitude_dir, (1.0, 0.0, 0.0))))
        assert angle_deg == pytest.approx(1.0)

    def test_a_sub_floor_demand_holds_fire_but_the_gimbal_keeps_tracking(self) -> None:
        from puffsat_sim.guidance import terminal_tick

        tick = terminal_tick(
            zem_est_m=(0.001, 0.0, 0.0),
            knowledge_sigma_m=0.0,
            t_go_s=100.0,
            feedforward_m_s2=(0.0, 0.0, 0.0),
            attitude_dir=(0.0, 1.0, 0.0),
            control_period_s=1.0,
            mass_kg=25.0,
        )

        assert not tick.fire
        angle_toward = np.degrees(np.arccos(np.dot(tick.attitude_dir, (0.0, 1.0, 0.0))))
        assert angle_toward == pytest.approx(1.0)

    def test_a_gated_out_estimate_is_engine_off_and_leaves_the_gimbal_parked(self) -> None:
        from puffsat_sim.guidance import terminal_tick

        tick = terminal_tick(
            zem_est_m=(10.0, 0.0, 0.0),
            knowledge_sigma_m=30.0,
            t_go_s=200.0,
            feedforward_m_s2=(0.0, 0.0, 0.0),
            attitude_dir=(0.0, 1.0, 0.0),
            control_period_s=1.0,
            mass_kg=25.0,
        )

        assert not tick.fire
        assert tick.thrust_n == 0.0
        assert tick.attitude_dir == (0.0, 1.0, 0.0)

    def test_the_drag_feedforward_shares_the_thruster_with_the_zem_command(self) -> None:
        from puffsat_sim.guidance import terminal_tick

        tick = terminal_tick(
            zem_est_m=(0.0, 0.0, 0.0),
            knowledge_sigma_m=0.0,
            t_go_s=100.0,
            feedforward_m_s2=(8e-4, 0.0, 0.0),
            attitude_dir=(1.0, 0.0, 0.0),
            control_period_s=1.0,
            mass_kg=25.0,
        )

        assert tick.fire
        assert tick.thrust_n == pytest.approx(0.02)
        assert np.allclose(tick.attitude_dir, (1.0, 0.0, 0.0))


def _fly_double_integrator(
    zem0: np.ndarray, span_s: float, dt_s: float, a_max_m_s2: float, gain: float = 3.0
) -> tuple[np.ndarray, dict[float, float]]:
    """ZOH ZEM loop on x'' = a; returns (final miss vector, {t_go: |ZEM|} history)."""
    from puffsat_sim.guidance import thrust_command, zem_acceleration

    mass = 25.0
    target = zem0.copy()
    x = np.zeros(3)
    v = np.zeros(3)
    history: dict[float, float] = {}
    steps = int(span_s / dt_s)
    for k in range(steps):
        t_go = span_s - k * dt_s
        zem = target - (x + v * t_go)
        history[t_go] = float(np.linalg.norm(zem))
        thrust_n, direction, _ = thrust_command(
            zem_acceleration((zem[0], zem[1], zem[2]), t_go, gain),
            mass_kg=mass,
            max_thrust_n=a_max_m_s2 * mass,
        )
        a = thrust_n / mass * np.asarray(direction)
        x = x + v * dt_s + 0.5 * a * dt_s**2
        v = v + a * dt_s
    return x - target, history


class TestZemLoopClosedForms:
    def test_unsaturated_loop_nulls_the_miss_and_zem_decays_as_t_go_cubed(self) -> None:
        zem0 = np.array([400.0, 0.0, 0.0])
        miss, history = _fly_double_integrator(zem0, span_s=180.0, dt_s=1.0, a_max_m_s2=1.0)

        assert np.linalg.norm(miss) < 0.01
        assert history[90.0] == pytest.approx(400.0 * (90.0 / 180.0) ** 3, rel=0.1)

    def test_saturated_loop_leaves_the_thrust_limited_residual(self) -> None:
        a_max = 0.016
        span = 180.0
        zem0 = np.array([4000.0, 0.0, 0.0])
        miss, _ = _fly_double_integrator(zem0, span_s=span, dt_s=1.0, a_max_m_s2=a_max)

        authority = 0.5 * a_max * span**2
        assert np.linalg.norm(miss) == pytest.approx(4000.0 - authority, rel=0.01)

    def test_with_tracker_noise_the_law_lands_within_the_measured_floor_multiple(self) -> None:
        from puffsat_sim.guidance import (
            CAPTURE_SIGMA_MAX_M,
            NavNoiseProcess,
            TrackerGrade,
            homing_floor_m,
            terminal_tick,
        )

        speed = 10_780.0
        a_max = 0.016
        span, dt, mass = 247.0, 1.0, 25.0

        def rms_lateral(sigma_theta: float, entry_m: float) -> float:
            lateral_sq = []
            for seed in range(8):
                rng = np.random.default_rng((20260612, seed))
                noise = NavNoiseProcess(TrackerGrade(sigma_theta, 1.0), rng)
                x = np.array([0.0, entry_m, 0.0])
                v = np.zeros(3)
                attitude: tuple[float, float, float] | None = None
                for k in range(int(span / dt)):
                    t_go = span - k * dt
                    range_m = speed * t_go
                    err = np.asarray(noise.sample((range_m, 0.0, 0.0), dt))
                    zem_est = -(x + err + v * t_go)
                    tick = terminal_tick(
                        (zem_est[0], zem_est[1], zem_est[2]),
                        knowledge_sigma_m=sigma_theta * range_m,
                        t_go_s=t_go,
                        feedforward_m_s2=(0.0, 0.0, 0.0),
                        attitude_dir=attitude,
                        control_period_s=dt,
                        mass_kg=mass,
                    )
                    attitude = tick.attitude_dir
                    a = tick.thrust_n / mass * np.asarray(attitude) if tick.fire else np.zeros(3)
                    x = x + v * dt + 0.5 * a * dt**2
                    v = v + a * dt
                lateral_sq.append(float(x[1] ** 2 + x[2] ** 2))
            return float(np.sqrt(np.mean(lateral_sq)))

        # The measured C3b structure: ~3x the closed-form floor at the 10 µrad design
        # point (the price of correlated knowledge error on a slew-limited engine);
        # capture-grade at 2 µrad.  The ungated law measured ~100x the floor.
        floor_10 = homing_floor_m(10e-6, speed, a_max)
        assert rms_lateral(10e-6, entry_m=0.0) < 5.0 * floor_10
        assert rms_lateral(10e-6, entry_m=400.0) < 5.0 * floor_10
        assert rms_lateral(2e-6, entry_m=0.0) < CAPTURE_SIGMA_MAX_M


class TestTrackerNoise:
    def test_angle_grade_scales_lateral_noise_with_range_and_range_noise_stays_fixed(
        self,
    ) -> None:
        from puffsat_sim.guidance import TrackerGrade, position_noise

        grade = TrackerGrade(sigma_theta_rad=10e-6, sigma_range_m=1.0)
        rng = np.random.default_rng(20260612)
        los = (700e3, 0.0, 0.0)
        draws = np.array([position_noise(grade, rng, los) for _ in range(4000)])

        assert float(np.std(draws[:, 0])) == pytest.approx(1.0, rel=0.05)
        assert float(np.std(draws[:, 1])) == pytest.approx(10e-6 * 700e3, rel=0.05)
        assert float(np.std(draws[:, 2])) == pytest.approx(10e-6 * 700e3, rel=0.05)

    def test_angle_grade_axes_follow_an_oblique_line_of_sight(self) -> None:
        from puffsat_sim.guidance import TrackerGrade, position_noise

        grade = TrackerGrade(sigma_theta_rad=50e-6, sigma_range_m=1.0)
        rng = np.random.default_rng(7)
        los = np.array([300e3, -400e3, 0.0])
        unit = los / np.linalg.norm(los)
        draws = np.array(
            [position_noise(grade, rng, (los[0], los[1], los[2])) for _ in range(4000)]
        )

        along = draws @ unit
        lateral = draws - np.outer(along, unit)
        sigma_lat = 50e-6 * 500e3
        assert float(np.std(along)) == pytest.approx(1.0, rel=0.05)
        assert float(np.mean(np.linalg.norm(lateral, axis=1))) == pytest.approx(
            sigma_lat * np.sqrt(np.pi / 2.0), rel=0.05
        )

    def test_constant_grade_is_isotropic_and_range_independent(self) -> None:
        from puffsat_sim.guidance import TrackerGrade, position_noise

        grade = TrackerGrade(sigma_theta_rad=None, sigma_range_m=1.0)
        rng = np.random.default_rng(11)
        draws = np.array([position_noise(grade, rng, (5_000e3, 0.0, 0.0)) for _ in range(4000)])

        assert np.allclose(np.std(draws, axis=0), 1.0, rtol=0.05)


class TestNavNoiseProcess:
    def test_knowledge_error_is_stationary_at_the_sigma_theta_r_envelope(self) -> None:
        from puffsat_sim.guidance import NavNoiseProcess, TrackerGrade

        grade = TrackerGrade(sigma_theta_rad=10e-6, sigma_range_m=1.0)
        process = NavNoiseProcess(grade, np.random.default_rng(3), tau_s=10.0)
        los = (700e3, 0.0, 0.0)
        draws = np.array([process.sample(los, dt_s=1.0) for _ in range(20000)])

        assert float(np.std(draws[:, 0])) == pytest.approx(1.0, rel=0.1)
        assert float(np.std(draws[:, 1])) == pytest.approx(10e-6 * 700e3, rel=0.1)
        assert float(np.std(draws[:, 2])) == pytest.approx(10e-6 * 700e3, rel=0.1)

    def test_consecutive_errors_are_correlated_over_the_smoothing_time(self) -> None:
        from puffsat_sim.guidance import NavNoiseProcess, TrackerGrade

        grade = TrackerGrade(sigma_theta_rad=10e-6, sigma_range_m=1.0)
        process = NavNoiseProcess(grade, np.random.default_rng(4), tau_s=10.0)
        los = (700e3, 0.0, 0.0)
        draws = np.array([process.sample(los, dt_s=1.0) for _ in range(20000)])

        lateral = draws[:, 1]
        lag1 = float(np.corrcoef(lateral[:-1], lateral[1:])[0, 1])
        assert lag1 == pytest.approx(np.exp(-1.0 / 10.0), abs=0.05)

    def test_the_envelope_shrinks_with_range(self) -> None:
        from puffsat_sim.guidance import NavNoiseProcess, TrackerGrade

        grade = TrackerGrade(sigma_theta_rad=10e-6, sigma_range_m=0.0)
        process = NavNoiseProcess(grade, np.random.default_rng(5), tau_s=10.0)
        far = np.array([process.sample((700e3, 0.0, 0.0), dt_s=1.0) for _ in range(5000)])
        near = np.array([process.sample((50e3, 0.0, 0.0), dt_s=1.0) for _ in range(5000)])

        assert float(np.std(near[:, 1])) == pytest.approx(
            float(np.std(far[:, 1])) * 50.0 / 700.0, rel=0.15
        )


class TestPredictedZem:
    def test_zem_is_zero_when_the_target_sits_on_the_ballistic_path(self) -> None:
        from puffsat_sim.estimation import two_body_j2_flow
        from puffsat_sim.guidance import predicted_zem

        state = np.array([6_378e3 + 800e3, 0.0, 0.0, 0.0, 7.5e3, 6.0e3])
        t_go = 170.0
        target = two_body_j2_flow(state, t_go)[:3]

        zem = predicted_zem(state, (target[0], target[1], target[2]), t_go)

        assert np.linalg.norm(zem) < 1e-6

    def test_a_displaced_target_reads_back_as_the_displacement(self) -> None:
        from puffsat_sim.estimation import two_body_j2_flow
        from puffsat_sim.guidance import predicted_zem

        state = np.array([6_378e3 + 800e3, 0.0, 0.0, 0.0, 7.5e3, 6.0e3])
        t_go = 170.0
        offset = np.array([120.0, -50.0, 80.0])
        target = two_body_j2_flow(state, t_go)[:3] + offset

        zem = predicted_zem(state, (target[0], target[1], target[2]), t_go)

        assert np.allclose(zem, offset, atol=1e-6)


class TestThrustCommand:
    def test_within_the_cap_thrust_is_mass_times_accel_along_the_command(self) -> None:
        thrust_n, direction, saturated = thrust_command(
            (3e-4, 0.0, 4e-4), mass_kg=25.0, max_thrust_n=0.4
        )

        assert thrust_n == 25.0 * 5e-4
        assert np.allclose(direction, (0.6, 0.0, 0.8))
        assert not saturated

    def test_caps_at_the_actuator_and_flags_saturation(self) -> None:
        thrust_n, direction, saturated = thrust_command(
            (0.03, 0.0, 0.04), mass_kg=25.0, max_thrust_n=0.4
        )

        assert thrust_n == 0.4
        assert np.allclose(direction, (0.6, 0.0, 0.8))
        assert saturated

    def test_zero_command_is_engine_off_not_a_nan_direction(self) -> None:
        thrust_n, direction, saturated = thrust_command(
            (0.0, 0.0, 0.0), mass_kg=25.0, max_thrust_n=0.4
        )

        assert thrust_n == 0.0
        assert np.all(np.isfinite(direction))
        assert not saturated

    def test_commands_below_the_proportional_floor_are_engine_off(self) -> None:
        thrust_n, _, saturated = thrust_command(
            (1e-4, 0.0, 0.0), mass_kg=25.0, max_thrust_n=0.4, floor_n=0.005
        )

        assert thrust_n == 0.0
        assert not saturated

    def test_commands_at_or_above_the_floor_pass_through(self) -> None:
        thrust_n, _, _ = thrust_command(
            (4e-4, 0.0, 0.0), mass_kg=25.0, max_thrust_n=0.4, floor_n=0.005
        )

        assert thrust_n == pytest.approx(0.01)
