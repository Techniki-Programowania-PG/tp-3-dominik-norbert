import scikit_build_example as sc
import os

def test_signal(signal_name, generator_fn, freq, duration, sampling_rate, cutoff):
    print(f"\n--- {signal_name} ---")
    print(f"Generowanie sygnalu {signal_name}...")
    signal = generator_fn(freq, duration, sampling_rate)

    num_samples = len(signal)
    t = [i / sampling_rate for i in range(num_samples)]

    print("Obliczanie DFT...")
    dft_result = sc.dft(signal)

    print("Wizualizacja sygnalu i DFT...")
    filename_dft = f"{signal_name}_dft.png"
    sc.visualize_signal_and_dft(signal, t, dft_result, filename_dft)
    print(f"Zapisano: {filename_dft}")

    print("Obliczanie IDFT...")
    reconstructed = sc.idft(dft_result)

    print("Wizualizacja DFT -> IDFT...")
    filename_idft = f"{signal_name}_idft.png"
    sc.visualize_dft_idft(signal, dft_result, reconstructed, filename_idft)
    print(f"Zapisano: {filename_idft}")

    print("Test rekonstrukcji...")
    reconstructed_real = [x.real for x in reconstructed]
    error = sum((s - r)**2 for s, r in zip(signal, reconstructed_real)) / len(signal)
    print(f"Blad sredniokwadratowy rekonstrukcji: {error:.6e}")

    filtered_dft = sc.low_pass_filter_1d(dft_result, sampling_rate, cutoff)
    filtered_signal = sc.idft(filtered_dft)

    filename_filtered = f"{signal_name}_filtered.png"
    sc.visualize_filtered_signal_1d(signal, filtered_dft, filtered_signal, filename_filtered)
    print(f"Zapisano: {filename_filtered}")


def main():
    print("Wprowadz parametry sygnalow:")
    freq = float(input("Czestotliwosc [Hz]: ") or 25)
    duration = float(input("Czas trwania [s]: ") or 0.5)
    sampling_rate = float(input("Czestotliwosc prbkowania [Hz]: ") or 2000) 
    cutoff = float(input(f"Podaj czestotliwosc odciecia (filtracji 1d) dla sygnalow  [Hz]: ") or 5)

    os.makedirs("output", exist_ok=True)
    os.chdir("output")

    test_signal("sinusoidalny", sc.generate_sin, freq, duration, sampling_rate, cutoff)
    test_signal("cosinusoidalny", sc.generate_cos, freq, duration, sampling_rate, cutoff)
    test_signal("prostokatny", sc.generate_square, freq, duration, sampling_rate, cutoff)
    test_signal("piloksztaltny", sc.generate_sawtooth, freq, duration, sampling_rate, cutoff)

    print("\nWszystkie sygnaly zostaly przetworzone. Zdjecia wynikow zostaly umieszczone w pliku tests/output.")


if __name__ == "__main__":
    main()
