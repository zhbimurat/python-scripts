# VIM for the critical load of a self-weighted column / ВИМ для критической силы стержня с собственным весом / Меншікті салмағы бар сырықтың сыни күші үшін ВИӘ

Supplementary material S1 to the article
*"Determining the critical buckling load of a self-weighted column using the variational iteration method"*
(Б. Ж. Сағындықов, Ж. Бимұрат).

**Languages:** [English](#english) · [Русский](#русский) · [Қазақша](#қазақша)

---

## English

### What the script does

`Supplementary_S1_vim_column.py` reproduces all numerical results of the article:

1. **Self-weight cantilever column** (clamped at the base, free at the top):
  `θ'' + (λ + μ(1−ξ))θ = 0, θ(0) = 0, θ'(1) = 0`, with `λ = P·l²/EI`, `μ = q·l³/EI`.
  - SymPy derives the first and second variational iterations (VIM-1, VIM-2) and their eigenvalue equations.
  - Eigenvalues are compared with a **reference high-accuracy numerical solution** (shooting method).
  - Greenhill's parameter `μ*` (reference ≈ 7.8373; VIM-1 → 6; VIM-2 → 8.2918) is computed.
  - Prints **Table 1** (μ = 0…7).
2. **Column on a Winkler foundation** (fourth-order equation):
  `Y'''' + (λ + μ(1−ξ))Y'' − μY' + κY = 0`, with `κ = k_w·l⁴/EI`.
  - SymPy builds the VIM-1 determinant `a11·a22 − a12·a21`.
  - Prints **Table 3** (VIM-1 vs. reference solution).
3. **Numerical example** (steel mine prop): prints λ, μ, κ, the critical forces (reference and VIM-1) and the classical Euler value.



### Requirements

- Python 3.9 or newer (tested with Python 3.12.3)
- Libraries from `requirements.txt` (tested with numpy 2.5.4, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0)



### Installation

```bash
git clone <repository-url>
cd <repository-folder>
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```



### Run

```bash
python Supplementary_S1_vim_column.py
```

The script has **no command-line arguments**. Runtime is about 1–2 minutes (the symbolic integration for the Winkler problem takes most of the time). Output goes to the console.

### Settings (edit in the source code)


| What                                | Where                             | Default                                                                                |
| ----------------------------------- | --------------------------------- | -------------------------------------------------------------------------------------- |
| mpmath precision for VIM roots      | `mp.mp.dps` (top of file)         | `40` digits                                                                            |
| ODE tolerances, self-weight problem | `shoot_selfweight()`              | `rtol=1e-12`, `atol=1e-12`, `DOP853`                                                   |
| ODE tolerances, Winkler problem     | `shoot_winkler_det()`             | `rtol=1e-10`, `atol=1e-12`                                                             |
| Root search window for VIM          | `first_root_mp(F1, m, 0.002, 10)` | λ ∈ (0.002; 10), excludes the second eigenvalue                                        |
| Values of μ, κ in the tables        | `__main__` block                  | μ = 0…7; Winkler list of (μ, κ) pairs                                                  |
| Input data of the example           | end of `__main__`                 | E = 2.1·10¹¹ Pa, l = 3 m, 0.1×0.1 m section, ρ = 7850 kg/m³, P = 350 kN, k_s = 6 MN/m³ |




### Expected output (key values)

- Greenhill parameter: reference `μ* = 7.8373`; VIM-1 `6`; VIM-2 `8.2918`.
- Table 1, row μ = 1: reference `2.1679`, VIM-1 `2.1093` (−2.70 %), VIM-2 `2.1687` (+0.04 %).
- Table 1, rows μ = 6 and 7: VIM-1 has no positive first root and is printed as `(<=0)`.
- Example: λ = 1.800, μ = 0.01188, κ = 27.77; reference P_cr ≈ 1289.0 kN, VIM-1 ≈ 1491.0 kN, Euler ≈ 479.8 kN.



### Notes and limitations

- The window λ < 10 prevents the second eigenvalue from being mistaken for the first one.
- For μ > μ* the first eigenvalue is negative (the column buckles under its own weight); VIM-1/VIM-2 are not reliable in that region and it is not part of this script.
- Winkler problem: only the first VIM iteration is implemented.
- The script is for method validation, not for design calculations.

## Русский

### Что делает скрипт

`Supplementary_S1_vim_column.py` воспроизводит все численные результаты статьи:

1. **Консольный стержень с собственным весом** (защемлён внизу, свободен вверху):
  `θ'' + (λ + μ(1−ξ))θ = 0, θ(0) = 0, θ'(1) = 0`, где `λ = P·l²/EI`, `μ = q·l³/EI`.
  - SymPy выводит первую и вторую итерации ВИМ (ВИМ-1, ВИМ-2) и уравнения для собственных значений.
  - Собственные значения сравниваются с **эталонным высокоточным численным решением** (метод пристрелки).
  - Вычисляется параметр Гринхилла `μ*` (эталон ≈ 7,8373; ВИМ-1 → 6; ВИМ-2 → 8,2918).
  - Выводится **таблица 1** (μ = 0…7).
2. **Стержень на основании Винклера** (уравнение четвёртого порядка):
  `Y'''' + (λ + μ(1−ξ))Y'' − μY' + κY = 0`, где `κ = k_w·l⁴/EI`.
  - SymPy строит определитель ВИМ-1 `a11·a22 − a12·a21`.
  - Выводится **таблица 3** (ВИМ-1 и эталонное решение).
3. **Числовой пример** (стальная шахтная стойка): выводятся λ, μ, κ, критические силы (эталонная и по ВИМ-1) и классическое значение Эйлера.



### Требования

- Python 3.9 и новее (проверено на Python 3.12.3)
- Библиотеки из `requirements.txt` (проверено на numpy 2.5.4, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0)



### Установка

```bash
git clone <ссылка-на-репозиторий>
cd <папка-репозитория>
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```



### Запуск

```bash
python Supplementary_S1_vim_column.py
```

**Аргументов командной строки нет.** Время работы около 1–2 минут (основное время занимает символьное интегрирование для задачи Винклера). Результат выводится в консоль.

### Настройки (меняются в коде)


| Что                                    | Где                               | По умолчанию                                                                           |
| -------------------------------------- | --------------------------------- | -------------------------------------------------------------------------------------- |
| Точность mpmath для корней ВИМ         | `mp.mp.dps` (начало файла)        | `40` знаков                                                                            |
| Допуски ОДУ, задача о собственном весе | `shoot_selfweight()`              | `rtol=1e-12`, `atol=1e-12`, `DOP853`                                                   |
| Допуски ОДУ, задача Винклера           | `shoot_winkler_det()`             | `rtol=1e-10`, `atol=1e-12`                                                             |
| Окно поиска корней ВИМ                 | `first_root_mp(F1, m, 0.002, 10)` | λ ∈ (0,002; 10), исключает второе собственное значение                                 |
| Значения μ, κ в таблицах               | блок `__main__`                   | μ = 0…7; список пар (μ, κ) для Винклера                                                |
| Исходные данные примера                | конец `__main__`                  | E = 2,1·10¹¹ Па, l = 3 м, сечение 0,1×0,1 м, ρ = 7850 кг/м³, P = 350 кН, k_s = 6 МН/м³ |




### Ожидаемый результат (ключевые значения)

- Параметр Гринхилла: эталон `μ* = 7,8373`; ВИМ-1 `6`; ВИМ-2 `8,2918`.
- Таблица 1, строка μ = 1: эталон `2,1679`, ВИМ-1 `2,1093` (−2,70 %), ВИМ-2 `2,1687` (+0,04 %).
- Таблица 1, строки μ = 6 и 7: у ВИМ-1 нет положительного первого корня, выводится `(<=0)`.
- Пример: λ = 1,800, μ = 0,01188, κ = 27,77; эталонная P_cr ≈ 1289,0 кН, ВИМ-1 ≈ 1491,0 кН, Эйлер ≈ 479,8 кН.



### Примечания и ограничения

- Окно λ < 10 не даёт принять второе собственное значение за первое.
- При μ > μ* первое собственное значение отрицательно (стержень теряет устойчивость под собственным весом); ВИМ-1/ВИМ-2 в этой области ненадёжны, и она в скрипт не входит.
- Задача Винклера: реализована только первая итерация ВИМ.
- Скрипт предназначен для проверки метода, а не для проектных расчётов.

## Қазақша



### Скрипттің міндеті

`Supplementary_S1_vim_column.py` мақаладағы барлық сандық нәтижелерді қайта өндіреді:

1. **Өз салмағы бар консольді сырық** (төменгі ұшы қатаң бекітілген, жоғарғы ұшы бос):
  `θ'' + (λ + μ(1−ξ))θ = 0, θ(0) = 0, θ'(1) = 0`; мұнда `λ = P·l²/EI`, `μ = q·l³/EI`.
  - SymPy вариациялық итерацияның бірінші және екінші қадамын (ВИӘ-1, ВИӘ-2) және меншікті мәндер теңдеулерін қорытады.
  - Меншікті мәндер **эталондық жоғары дәлдіктегі сандық шешіммен** (соққы әдісі) салыстырылады.
  - Гринхилл параметрі `μ*` есептеледі (эталон ≈ 7,8373; ВИӘ-1 → 6; ВИӘ-2 → 8,2918).
  - **1-кесте** шығарылады (μ = 0…7).
2. **Винклер іргетасындағы сырық** (төртінші ретті теңдеу):
  `Y'''' + (λ + μ(1−ξ))Y'' − μY' + κY = 0`; мұнда `κ = k_w·l⁴/EI`.
  - SymPy ВИӘ-1 анықтауышын `a11·a22 − a12·a21` құрады.
  - **3-кесте** шығарылады (ВИӘ-1 және эталондық шешім).
3. **Сандық мысал** (болат шахта тіреуіші): λ, μ, κ, сыни күштер (эталондық және ВИӘ-1 бойынша) және классикалық Эйлер мәні шығарылады.



### Талаптар

- Python 3.9 және одан жоғары (Python 3.12.3 нұсқасында тексерілген)
- `requirements.txt` ішіндегі кітапханалар (numpy 2.5.4, scipy 1.18.1, sympy 1.14.0, mpmath 1.3.0 нұсқаларында тексерілген)



### Орнату

```bash
git clone <репозиторий-сілтемесі>
cd <репозиторий-қалтасы>
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```



### Іске қосу

```bash
python Supplementary_S1_vim_column.py
```

**Пәрмен жолы аргументтері жоқ.** Орындалу уақыты шамамен 1–2 минут (негізгі уақытты Винклер есебі үшін символдық интегралдау алады). Нәтиже консольге шығады.

### Баптаулар (кодта өзгертіледі)


| Не                                      | Қайда                             | Әдепкі                                                                                |
| --------------------------------------- | --------------------------------- | ------------------------------------------------------------------------------------- |
| ВИӘ түбірлеріне арналған mpmath дәлдігі | `mp.mp.dps` (файл басы)           | `40` таңба                                                                            |
| ҚДТ рұқсат шамалары, өз салмағы есебі   | `shoot_selfweight()`              | `rtol=1e-12`, `atol=1e-12`, `DOP853`                                                  |
| ҚДТ рұқсат шамалары, Винклер есебі      | `shoot_winkler_det()`             | `rtol=1e-10`, `atol=1e-12`                                                            |
| ВИӘ түбірін іздеу терезесі              | `first_root_mp(F1, m, 0.002, 10)` | λ ∈ (0,002; 10), екінші меншікті мәнді алып тастайды                                  |
| Кестелердегі μ, κ мәндері               | `__main__` блогы                  | μ = 0…7; Винклер үшін (μ, κ) жұптарының тізімі                                        |
| Мысалдың бастапқы деректері             | `__main__` соңы                   | E = 2,1·10¹¹ Па, l = 3 м, қимасы 0,1×0,1 м, ρ = 7850 кг/м³, P = 350 кН, k_s = 6 МН/м³ |




### Күтілетін нәтиже (негізгі мәндер)

- Гринхилл параметрі: эталон `μ* = 7,8373`; ВИӘ-1 `6`; ВИӘ-2 `8,2918`.
- 1-кесте, μ = 1 жолы: эталон `2,1679`, ВИӘ-1 `2,1093` (−2,70 %), ВИӘ-2 `2,1687` (+0,04 %).
- 1-кесте, μ = 6 және 7 жолдары: ВИӘ-1 үшін оң бірінші түбір жоқ, `(<=0)` шығарылады.
- Мысал: λ = 1,800, μ = 0,01188, κ = 27,77; эталондық P_cr ≈ 1289,0 кН, ВИӘ-1 ≈ 1491,0 кН, Эйлер ≈ 479,8 кН.



### Ескертпелер және шектеулер

- λ < 10 терезесі екінші меншікті мәнді біріншісі деп қабылдауға жол бермейді.
- μ > μ* болғанда бірінші меншікті мән теріс (сырық өз салмағынан орнықтылығын жоғалтады); бұл аймақта ВИӘ-1/ВИӘ-2 сенімсіз, сондықтан ол скриптке кірмейді.
- Винклер есебі: ВИӘ-нің тек бірінші итерациясы жүзеге асырылған.
- Скрипт әдісті тексеруге арналған, жобалық есептеулерге емес.



