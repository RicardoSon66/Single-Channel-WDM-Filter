# **Ring Resonator 기반 단일채널 WDM Filter 구현 및 시뮬레이션**  
이는 이전 프로젝트에서 문헌 근사치(ng≈4.2)로 계산했던 18.2nm와 거의 비슷한 값이지만, 이번에는 근사가 아니라 실제 스펙(radius=5.0, width=0.4, thickness=0.22)에 대한 모드솔빙 결과이므로 더 신뢰할 수 있는 예측값입니다.

## 2-2. Add-Drop 스펙트럼 형태 예측 (TMM)

이전 프로젝트에서 사용했던 all-pass 링 공식

$$T = \frac{a^2 - 2ra\cos\theta + r^2}{1 - 2ra\cos\theta + (ra)^2}$$

은 through 포트 하나만 있는 all-pass 구조에 대한 식이었습니다. 이번에는 drop 포트가 추가된 add-drop 구조이므로, 버스-링 결합계수를 through측 t1, drop측 t2로 분리한 아래 형태를 사용합니다.

$$T_{through}=\left|\frac{t_1-t_2 a e^{j\phi}}{1-t_1 t_2 a e^{j\phi}}\right|^2,\qquad T_{drop}=\left|\frac{-\kappa_1\kappa_2\sqrt{a}\,e^{j\phi/2}}{1-t_1 t_2 a e^{j\phi}}\right|^2 \quad(\kappa_i=\sqrt{1-t_i^2})$$

t1, t2(결합), a(loss)는 gap=0.15μm만으로 해석적으로 바로 구하기는 어려우므로 별도의 directional coupler 모드해석이 필요합니다. 따라서 이 값들은 3장의 FDTD 결과를 얻은 뒤, 4장에서 역산하여 채워 넣었습니다. 즉 이번 프로젝트의 이론적 예측은 FSR은 시뮬레이션 전에 예측, 스펙트럼 형태는 FDTD 결과와 함께 검증하는 2단계로 진행했습니다.

# **3. 시뮬레이션 세팅**

Source는 `ModeSource`를 사용했고 λ0=1.55μm 중심으로 wavelength_width=0.1μm(1.50~1.60μm)의 GaussianPulse입니다. Through와 Drop 두 포트에 각각 `FluxMonitor`를 배치했고, 링 내부에 `FieldTimeMonitor`를 두어 정상상태 필드도 함께 관찰했습니다.

이전 Meep 기반 시뮬레이션과 다르게, 이번에는 별도의 reference straight waveguide의 정규화 과정 없이 Tidy3D의 flux 결과를 바로 dB로 환산해서 사용했습니다(`Decibels = 10·log10(P/P_in)`). grid는 `auto(wavelength=1.55, min_steps_per_wvl=15)`, `run_time=30ps`, 경계조건은 전면 PML입니다.

```python
modesource = td.ModeSource(
    center=(-5, -wg_center_y, 0), size=(0, 1.6, 1.32),
    source_time=td.GaussianPulse(freq0=freq0, fwidth=fwidth),
    direction='+', mode_spec=mode_spec, mode_index=0,
)
ThroughMonitor = td.FluxMonitor(center=(5, -wg_center_y, 0), size=(0, 1.6, 1.32), freqs=freqs, name='ThroughMonitor')
DropMonitor    = td.FluxMonitor(center=(-5, wg_center_y, 0), size=(0, 1.6, 1.32), freqs=freqs, name='DropMonitor')
```


# **4. 스펙트럼 분석 및 동작 파장 선정**

Through/Drop 스펙트럼을 추출한 결과는 다음과 같습니다.

![FDTD Through 스펙트럼](./ThroughMonitor.png)
![FDTD Drop 스펙트럼](./DropMonitor.png)

이전 원형 Ring Resonator(gap=0.2μm)에서는 through의 dip이 -0.02~-0.1dB 수준으로 거의 관측이 불가능한 수준이었는데, 이번 스펙(gap=0.15μm, 3D FDTD)에서는 FSR 간격으로 뚜렷한 dip이 여러 개 관측되었습니다.  
이 스펙트럼에서 필터 동작을 확인할 두 파장, 즉 drop 채널로 쓰는 공진 파장 λ1과 비교 기준이 되는 비공진 파장 λ2를 다음과 같이 선정했습니다.

| 구분 | 파장 | Through | Drop |
|---|---|---|---|
| λ1 (공진, drop 채널) | 1.56774 μm | −20.69708 dB | −0.2974766 dB |
| λ2 (비공진, 기준 파장) | 1.54048 μm | −0.0234321 dB | −24.55664 dB |

λ1은 through에서 깊게 꺼지고 drop으로 대부분 에너지가 전달되는 공진 파장, λ2는 반대로 through를 그대로 통과하고 drop으로는 거의 새어나가지 않는 비공진 파장입니다. 따라서 링 하나가 λ1 채널 하나를 골라내는 단일 채널 add-drop 필터로 동작함을 확인했습니다. 다만 링은 FSR마다 공진하므로 λ1에서 FSR 간격만큼 떨어진 파장도 함께 drop됩니다.

`ResonanceFinder`로 추출한 Q값은 다음과 같습니다.

| freq (Hz) | wavelength (μm) | Q |
|---|---|---|
| 1.867731e+14 | 1.605116 | 771.7 |
| 1.867738e+14 | 1.605109 | 1281.0 |
| 1.889927e+14 | 1.586265 | 1501.8 |
| **1.912390e+14** | **1.567633 (≈λ1)** | **1754.8** |
| 1.934702e+14 | 1.549554 | 2122.0 |
| 1.956994e+14 | 1.531903 | 2477.3 |
| 1.979503e+14 | 1.514484 | 2972.6 |

이전 원형 Ring Resonator의 Harminv 결과(Q=3112.8, 단일 모드)와 비교하면, 이번에는 대역 전체에서 여러 개의 모드가 규칙적으로 검출되었고 소스 대역(1.50~1.60 μm) 안의 공진에서 Q값은 약 1500~3000으로 합리적인 범위입니다(1.605 μm 부근 값은 소스 대역 밖이라 신뢰하기 어렵습니다). Q가 파장이 짧아질수록(고주파일수록) 커지는 경향도 보이는데, 이는 손실/결합 조건이 파장에 따라 서서히 변하기 때문으로 해석됩니다.

정상상태(t=3.22ps) 필드 분포도 확인했습니다. λ1, λ2를 동시에 여기했을 때 링 내부에 정상적인 공진 모드 형태가 형성됨을 직접 눈으로 확인할 수 있었습니다.

![정상상태 필드 분포](./animation.gif)

# **5. 예측값과 비교 검증**

## 5-1. FSR 검증

| | 값 |
|---|---|
| 2절 이론 예측 FSR | 18.071 nm |
| FDTD 실측 FSR (인접 공진 λ 간격, 1.567633→1.549554μm) | 18.079 nm |
| 오차 | **0.008 nm (≈0.04%)** |

Full 3D FDTD를 돌리기 전, 단면 모드솔빙만으로 예측한 FSR이 실제 FDTD 결과와 사실상 일치했습니다.

## 5-2. TMM 이론 스펙트럼 vs FDTD 실측 비교 및 비공진 파장(λ2) 교차 검증

2-2절의 t1, t2, a는 λ1(공진점) 한 지점의 through/drop/Q 값으로만 역산(least-squares fitting)했습니다.
