/**
 * @file DatabaseManager.h
 * @brief Centralized database management for OpenMagnetics components
 * 
 * This class encapsulates all database operations for cores, materials,
 * wires, bobbins, and insulation materials, replacing the global variables
 * with a proper singleton that supports testing and thread-safe access.
 */

#pragma once

#include <MAS.hpp>
#include <map>
#include <vector>
#include <string>
#include <mutex>
#include <optional>
#include <memory>

// Forward declarations to avoid circular dependencies
namespace OpenMagnetics {
    class Core;
    class Wire;
    class Bobbin;
    class InsulationMaterial;
}

namespace OpenMagnetics {

/**
 * @brief Centralized manager for component databases
 * 
 * This singleton class manages all component databases including cores,
 * materials, shapes, wires, bobbins, and insulation materials. It provides
 * thread-safe access and supports database reloading for testing.
 * 
 * @note Use DatabaseManager::getInstance() to access the singleton.
 * 
 * Example usage:
 * @code
 *   auto& db = DatabaseManager::getInstance();
 *   db.loadCores();
 *   auto core = db.findCoreByName("E55_21");
 * @endcode
 */
class DatabaseManager {
public:
    /**
     * @brief Get the singleton instance
     * @return Reference to the DatabaseManager instance
     */
    static DatabaseManager& getInstance();
    
    // Prevent copying and moving
    DatabaseManager(const DatabaseManager&) = delete;
    DatabaseManager& operator=(const DatabaseManager&) = delete;
    DatabaseManager(DatabaseManager&&) = delete;
    DatabaseManager& operator=(DatabaseManager&&) = delete;
    
    // ========================================================================
    // Database Loading
    // ========================================================================
    
    /**
     * @brief Load all databases with default settings
     * @param withAliases Whether to load shape aliases
     * @param addInternalData Whether to add internal computed data
     */
    void loadAll(bool withAliases = true, bool addInternalData = true);
    
    /**
     * @brief Load cores database from embedded or external file
     * @param filePath Optional path to external file (uses embedded if empty)
     */
    void loadCores(const std::optional<std::string>& filePath = std::nullopt);
    
    /**
     * @brief Load core materials database
     * @param filePath Optional path to external file
     */
    void loadCoreMaterials(const std::optional<std::string>& filePath = std::nullopt);
    
    /**
     * @brief Load advanced core materials from external file
     * @param filePath Path to external file
     * @param onlyDataFromManufacturer Filter to manufacturer data only
     */
    void loadAdvancedCoreMaterials(const std::string& filePath, bool onlyDataFromManufacturer = true);
    
    /**
     * @brief Load core shapes database
     * @param withAliases Whether to include shape aliases
     * @param filePath Optional path to external file
     */
    void loadCoreShapes(bool withAliases = true, const std::optional<std::string>& filePath = std::nullopt);
    
    /**
     * @brief Load wires database
     * @param filePath Optional path to external file
     */
    void loadWires(const std::optional<std::string>& filePath = std::nullopt);
    
    /**
     * @brief Load bobbins database
     */
    void loadBobbins();
    
    /**
     * @brief Load insulation materials database
     */
    void loadInsulationMaterials();
    
    /**
     * @brief Load wire materials database
     */
    void loadWireMaterials();
    
    // ========================================================================
    // Database Clearing
    // ========================================================================
    
    /**
     * @brief Clear all loaded databases
     */
    void clearAll();
    
    /**
     * @brief Clear cores database only
     */
    void clearCores();
    
    /**
     * @brief Reset the database manager (for testing)
     * @note This clears all data and resets the internal state
     */
    void reset();
    
    // ========================================================================
    // Core Accessors
    // ========================================================================
    
    /**
     * @brief Find a core by name
     * @param name Core name to search for
     * @return The found Core
     * @throws CoreNotFoundException if core is not found
     */
    Core findCoreByName(const std::string& name) const;
    
    /**
     * @brief Try to find a core by name
     * @param name Core name to search for
     * @return Optional containing the Core if found
     */
    std::optional<Core> tryFindCoreByName(const std::string& name) const;
    
    /**
     * @brief Get all loaded cores
     * @return Reference to the cores vector
     */
    const std::vector<Core>& getCores() const { return _cores; }
    
    /**
     * @brief Get mutable access to cores (for modification)
     * @return Reference to the cores vector
     */
    std::vector<Core>& getMutableCores() { return _cores; }
    
    // ========================================================================
    // Material Accessors
    // ========================================================================
    
    /**
     * @brief Find a core material by name
     * @param name Material name
     * @return The found CoreMaterial
     * @throws CoreMaterialNotFoundException if not found
     */
    MAS::CoreMaterial findCoreMaterialByName(const std::string& name) const;
    
    /**
     * @brief Try to find a core material by name
     * @param name Material name
     * @return Optional containing the material if found
     */
    std::optional<MAS::CoreMaterial> tryFindCoreMaterialByName(const std::string& name) const;
    
    /**
     * @brief Get all core material names
     * @param manufacturer Optional filter by manufacturer
     * @return Vector of material names
     */
    std::vector<std::string> getCoreMaterialNames(const std::optional<std::string>& manufacturer = std::nullopt) const;
    
    /**
     * @brief Get all core materials
     * @param manufacturer Optional filter by manufacturer
     * @return Vector of materials
     */
    std::vector<MAS::CoreMaterial> getCoreMaterials(const std::optional<std::string>& manufacturer = std::nullopt) const;
    
    // ========================================================================
    // Shape Accessors
    // ========================================================================
    
    /**
     * @brief Find a core shape by name
     * @param name Shape name
     * @return The found CoreShape
     * @throws CoreShapeNotFoundException if not found
     */
    MAS::CoreShape findCoreShapeByName(const std::string& name) const;
    
    /**
     * @brief Get all core shape names
     * @return Vector of shape names
     */
    std::vector<std::string> getCoreShapeNames() const;
    
    /**
     * @brief Get shape names for a specific family
     * @param family Shape family
     * @return Vector of shape names
     */
    std::vector<std::string> getCoreShapeNames(MAS::CoreShapeFamily family) const;
    
    /**
     * @brief Get all shape families present in database
     * @return Vector of shape families
     */
    std::vector<MAS::CoreShapeFamily> getCoreShapeFamilies() const;
    
    /**
     * @brief Get all core shapes
     * @param includeToroidal Whether to include toroidal shapes
     * @return Vector of shapes
     */
    std::vector<MAS::CoreShape> getCoreShapes(bool includeToroidal = true) const;
    
    // ========================================================================
    // Wire Accessors  
    // ========================================================================
    
    /**
     * @brief Find a wire by name
     * @param name Wire name
     * @return The found Wire
     * @throws WireNotFoundException if not found
     */
    Wire findWireByName(const std::string& name) const;
    
    /**
     * @brief Find a wire by dimension
     * @param dimension Wire dimension
     * @param wireType Optional wire type filter
     * @param wireStandard Optional wire standard filter
     * @param obfuscate Whether to obfuscate the result
     * @return The found Wire
     */
    Wire findWireByDimension(double dimension,
                             std::optional<MAS::WireType> wireType = std::nullopt,
                             std::optional<MAS::WireStandard> wireStandard = std::nullopt,
                             bool obfuscate = true) const;
    
    /**
     * @brief Get all wire names
     * @return Vector of wire names
     */
    std::vector<std::string> getWireNames() const;
    
    /**
     * @brief Get all wires
     * @param wireType Optional filter by wire type
     * @param wireStandard Optional filter by wire standard
     * @return Vector of wires
     */
    std::vector<Wire> getWires(std::optional<MAS::WireType> wireType = std::nullopt,
                               std::optional<MAS::WireStandard> wireStandard = std::nullopt) const;
    
    // ========================================================================
    // Bobbin Accessors
    // ========================================================================
    
    /**
     * @brief Find a bobbin by name
     * @param name Bobbin name
     * @return The found Bobbin
     */
    Bobbin findBobbinByName(const std::string& name) const;
    
    /**
     * @brief Get all bobbin names
     * @return Vector of bobbin names
     */
    std::vector<std::string> getBobbinNames() const;
    
    /**
     * @brief Get all bobbins
     * @return Vector of bobbins
     */
    std::vector<Bobbin> getBobbins() const;
    
    // ========================================================================
    // Insulation Material Accessors
    // ========================================================================
    
    /**
     * @brief Find an insulation material by name
     * @param name Material name
     * @return The found InsulationMaterial
     */
    InsulationMaterial findInsulationMaterialByName(const std::string& name) const;
    
    /**
     * @brief Get all insulation material names
     * @return Vector of material names
     */
    std::vector<std::string> getInsulationMaterialNames() const;
    
    /**
     * @brief Get all insulation materials
     * @return Vector of materials
     */
    std::vector<InsulationMaterial> getInsulationMaterials() const;
    
    // ========================================================================
    // Wire Material Accessors
    // ========================================================================
    
    /**
     * @brief Find a wire material by name
     * @param name Material name
     * @return The found WireMaterial
     */
    MAS::WireMaterial findWireMaterialByName(const std::string& name) const;
    
    /**
     * @brief Get all wire material names
     * @return Vector of material names
     */
    std::vector<std::string> getWireMaterialNames() const;
    
    /**
     * @brief Get all wire materials
     * @return Vector of materials
     */
    std::vector<MAS::WireMaterial> getWireMaterials() const;
    
    // ========================================================================
    // Status Checking
    // ========================================================================
    
    /**
     * @brief Check if cores database is loaded
     * @return true if cores are loaded
     */
    bool isCoresLoaded() const { return !_cores.empty(); }
    
    /**
     * @brief Check if materials database is loaded
     * @return true if materials are loaded
     */
    bool isMaterialsLoaded() const { return !_coreMaterials.empty(); }
    
    /**
     * @brief Check if shapes database is loaded
     * @return true if shapes are loaded
     */
    bool isShapesLoaded() const { return !_coreShapes.empty(); }
    
    /**
     * @brief Check if wires database is loaded
     * @return true if wires are loaded
     */
    bool isWiresLoaded() const { return !_wires.empty(); }
    
private:
    // Private constructor for singleton
    DatabaseManager() = default;
    ~DatabaseManager() = default;
    
    // Database storage
    std::vector<Core> _cores;
    std::map<std::string, MAS::CoreMaterial> _coreMaterials;
    std::map<std::string, MAS::CoreShape> _coreShapes;
    std::vector<MAS::CoreShapeFamily> _coreShapeFamiliesInDatabase;
    std::map<std::string, Wire> _wires;
    std::map<std::string, Bobbin> _bobbins;
    std::map<std::string, InsulationMaterial> _insulationMaterials;
    std::map<std::string, MAS::WireMaterial> _wireMaterials;
    
    // Thread safety
    mutable std::mutex _mutex;
    
    // Loading state flags
    bool _coresLoaded = false;
    bool _materialsLoaded = false;
    bool _shapesLoaded = false;
    bool _wiresLoaded = false;
    bool _bobbinsLoaded = false;
    bool _insulationMaterialsLoaded = false;
    bool _wireMaterialsLoaded = false;
};

/**
 * @brief RAII helper for database transaction/testing
 * 
 * This class saves the current database state and restores it
 * when destroyed. Useful for testing scenarios.
 */
class DatabaseTransaction {
public:
    DatabaseTransaction();
    ~DatabaseTransaction();
    
    // Prevent copying
    DatabaseTransaction(const DatabaseTransaction&) = delete;
    DatabaseTransaction& operator=(const DatabaseTransaction&) = delete;
    
private:
    // Saved state would be stored here
    bool _committed = false;
};

} // namespace OpenMagnetics
