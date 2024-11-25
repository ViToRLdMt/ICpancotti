def convert_fe_ga203(input_file, output_file):
    """
    Converte arquivos Fe on Ga203 aplicando subtrações iterativas com base no terceiro valor do bloco.
    Mantém o valor base da segunda coluna na primeira linha de cada bloco.
    Remove os espaços entre as linhas no arquivo de saída.

    Args:
        input_file (str): Caminho do arquivo de entrada (.txt).
        output_file (str): Caminho do arquivo de saída (.txt).
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        base_values = None  # Para armazenar os valores base do bloco
        current_value = None  # Valor corrente para subtrações progressivas
        is_first_line = True  # Flag para controlar a primeira linha do bloco

        for line in infile:
            parts = line.strip().split()

            # Mantém as linhas de texto ou cabeçalhos sem alteração
            if any(c.isalpha() for c in line):
                outfile.write(line)
                continue

            # Verifica se a linha contém os valores base (início de um bloco)
            if len(parts) >= 3 and all(p.lstrip('-').replace('.', '', 1).isdigit() for p in parts[:3]):
                # Define os valores base do bloco
                base_values = [float(parts[0]), float(parts[1]), float(parts[2])]
                current_value = base_values[0]  # Inicializa o valor corrente com o primeiro valor do bloco
                is_first_line = True  # Reseta o controle de primeira linha
                # Escreve a linha original do cabeçalho
                outfile.write(line)
                continue

            # Caso a linha contenha um único número válido (dados do bloco)
            if len(parts) == 1 and parts[0].replace('.', '', 1).isdigit():
                if base_values is not None and current_value is not None:
                    # Pega o valor original
                    original_value = int(float(parts[0]))

                    # Mantém o valor base na primeira linha
                    if is_first_line:
                        outfile.write(f"{original_value:<10}{current_value:<15.5f}\n")
                        is_first_line = False  # Desativa o controle de primeira linha
                    else:
                        # Faz a subtração corretamente a partir da segunda linha
                        current_value += base_values[2]
                        outfile.write(f"{original_value:<10}{current_value:<15.5f}\n")
                continue

            # Mantém as linhas inválidas ou vazias como estão
            outfile.write(line)

    print(f"Conversão concluída! Arquivo salvo em: {output_file}")


# Inputs mantidos
input_file = input("Digite o caminho do arquivo de entrada: ").strip()
output_file = input("Digite o caminho do arquivo de saída: ").strip()

convert_fe_ga203(input_file, output_file)
