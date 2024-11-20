/*---------------------------------------------------------------------------*/
// * This file contains a function for read a single value.
/*---------------------------------------------------------------------------*/


#ifndef _READVALUECUH
#define _READVALUECUH

#pragma once



/** Read data from external file */
real_t readValueFromFile(std::string nombreArchivo) {
    real_t valor = 0.0f; // Inicializamos la variable valor con 0.0

    // Abre el archivo en modo lectura
    std::ifstream archivo(nombreArchivo);

    // Verifica si el archivo se ha abierto correctamente
    if (!archivo.is_open()) {
        std::cerr << "Error al abrir el archivo." << std::endl;
        // Si no se puede abrir el archivo, devuelve 0.0
        return valor;
    }

    // Lee el valor desde el archivo
    archivo >> valor;

    // Cierra el archivo después de leer
    archivo.close();

    // Devuelve el valor leído
    return valor;
}

#endif // -> #ifdef _READVALUECUH