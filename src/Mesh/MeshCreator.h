#pragma once

// Project files
#include <MacrosAndTypedefs.h>
#include <Miscellaneous.h>
#include <Parameter.h>

// deal.II files
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

/**
 * MeshCreator is a struct with only functions to create mesh of required
 * geometry. This struct has a static function to create mesh based on the
 * passed type. It includes also functionality to setup boundary IDs.
 */
struct MeshCreator {
  using MeshCreatorFunction = std::function<void(DomainParallelTriaTypeDistributed &)>;

  static void createMesh(DomainParallelTriaTypeDistributed &tria, const std::string &type) {
    _getCreator(type)(tria);
  }

private:
  static MeshCreatorFunction _getCreator(const std::string &type) {
    static const std::map<std::string, MeshCreatorFunction> creators = {
        {"hyper_rectangle",
         [](DomainParallelTriaTypeDistributed &tria) { _CreateHyperRectangleDomain(tria); }},
        {"hyper_L", [](DomainParallelTriaTypeDistributed &tria) { _CreateHyperLDomain(tria); }},
        {"abaqus", [](DomainParallelTriaTypeDistributed &tria) { _CreateAbaqusDomain(tria); }}};

    auto it = creators.find(type);
    if (it != creators.end()) {
      return it->second;
    } else {
      throw std::invalid_argument("Unknown mesh type: " + type);
    }
  }

  static void _CreateHyperRectangleDomain(DomainParallelTriaTypeDistributed &tria);
  static void _CreateHyperLDomain(DomainParallelTriaTypeDistributed &tria);
  static void _CreateAbaqusDomain(DomainParallelTriaTypeDistributed &tria);

  // Padded domains to avoid boundary effects in the filter
  static void _RemovePaddedCellsForBoundaryConditions(DomainParallelTriaTypeDistributed &tria);
  static void _CreatePaddedHyperRectangleDomain(DomainParallelTriaTypeDistributed &tria);
  static void _CreatePaddedHyperLDomain(DomainParallelTriaTypeDistributed &tria);
};