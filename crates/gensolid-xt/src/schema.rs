//! Exact ABC/Onshape 30.1 embedded-schema adapter.
//!
//! Wire definitions come from the pinned MIT/Apache-2.0 `parasolid-core` 0.3.7
//! dependency. Semantic roles and the 13006 type-table bound are checked against
//! Siemens' public JT v10 Rev-B, Appendix E.2 and E.5--E.7:
//! <https://www.plm.automation.siemens.com/en_us/Images/JT-v10-file-format-reference-rev-B_tcm1023-233786.pdf>.
//! Only the exact key found in the public ABC sample is enabled. Embedded edits
//! must copy every semantic field from the reviewed base; input names never
//! grant a newly inserted field a topology or geometry role.

use std::collections::{BTreeMap, BTreeSet};

use gensolid_core::{Error, Result};
use parasolid_core::{
    FieldValue, RawNode, SchemaEdit, SchemaKey, SchemaProvider, SchemaResolution, SchemaSource,
    SchemaTypeLookup, TypeDefinition,
};

pub(crate) const ABC_KEY: &str = "SCH_3001168_30100_13006";
const BASE: &str = "13006";
// Siemens Appendix E.6 lists the 13006 node-type table through HELIX_SU_FORM (184).
const BASE_TYPE_BOUND: u16 = 184;

pub(crate) struct AbcSchema {
    definitions: BTreeMap<u16, TypeDefinition>,
    semantic_fields: BTreeMap<u16, BTreeSet<String>>,
}

impl AbcSchema {
    pub(crate) fn new() -> Result<Self> {
        let profile = parasolid_core::schema::profiles::onshape_sch13006()
            .map_err(super::import::parser_error)?;
        let mut definitions = BTreeMap::new();
        let mut semantic_fields = BTreeMap::new();
        for original in profile.definitions() {
            let mut definition = original.clone();
            let mut roles = BTreeSet::new();
            if let Some(name) = type_role(definition.node_type) {
                definition.name = name.to_owned();
            }
            for field in &mut definition.fields {
                if let Some(role) = field_role(definition.node_type, &field.name) {
                    field.name = role.to_owned();
                    roles.insert(field.name.clone());
                }
            }
            semantic_fields.insert(definition.node_type, roles);
            definitions.insert(definition.node_type, definition);
        }
        Ok(Self {
            definitions,
            semantic_fields,
        })
    }

    /// Prevent deleted, inserted or duplicate input fields from impersonating a base role.
    pub(crate) fn audit(&self, resolutions: &[SchemaResolution]) -> Result<()> {
        for resolution in resolutions {
            let kind = resolution.definition.node_type;
            let Some(base) = self.definitions.get(&kind) else {
                continue;
            };
            let roles = self
                .semantic_fields
                .get(&kind)
                .ok_or_else(|| Error::internal("missing XT semantic field set"))?;
            match resolution.definition.source {
                SchemaSource::EmbeddedUnchanged => {
                    if resolution.definition.fields != base.fields {
                        return Err(Error::parse(format!(
                            "XT type {kind} changed an unchanged base definition"
                        )));
                    }
                }
                SchemaSource::EmbeddedDelta => {
                    let mut cursor = 0;
                    let mut copied = BTreeSet::new();
                    for edit in &resolution.edits {
                        match edit {
                            SchemaEdit::Copy { .. } => {
                                let field = base.fields.get(cursor).ok_or_else(|| {
                                    Error::parse("XT schema copy exceeds base definition")
                                })?;
                                if roles.contains(&field.name) {
                                    copied.insert(field.name.clone());
                                }
                                cursor += 1;
                            }
                            SchemaEdit::Delete { .. } => cursor += 1,
                            SchemaEdit::Insert { field, .. } | SchemaEdit::Append { field, .. } => {
                                if roles.contains(&field.name) {
                                    return Err(Error::unimplemented(format!(
                                        "XT type {kind} inserts semantic field {}; its meaning is not established",
                                        field.name
                                    )));
                                }
                            }
                            SchemaEdit::End { .. } => {}
                        }
                    }
                    if &copied != roles {
                        return Err(Error::unimplemented(format!(
                            "XT type {kind} deletes or replaces a semantic base field"
                        )));
                    }
                }
                _ => {
                    return Err(Error::unimplemented(format!(
                        "XT type {kind} did not retain its reviewed embedded-base provenance"
                    )));
                }
            }
        }
        Ok(())
    }

    /// These body extensions describe geometry the face-only target cannot retain.
    /// The reviewed ABC parts transmit null pointers for them. A file with the
    /// same key but convergent or nested-body data must not become a partial solid.
    pub(crate) fn audit_nodes(&self, nodes: &[RawNode]) -> Result<()> {
        for node in nodes.iter().filter(|node| node.node_type == 12) {
            for field in &node.fields {
                if matches!(
                    field.definition.name.as_str(),
                    "mesh"
                        | "polyline"
                        | "boundary_mesh"
                        | "boundary_polyline"
                        | "child"
                        | "mesh_offset_data"
                ) && field
                    .values
                    .iter()
                    .any(|value| matches!(value, FieldValue::PointerIndex(index) if *index != 0))
                {
                    return Err(Error::unimplemented(format!(
                        "XT body node {} has unsupported non-null {} geometry",
                        node.index, field.definition.name
                    )));
                }
            }
        }
        Ok(())
    }
}

impl SchemaProvider for AbcSchema {
    fn contains_schema(&self, schema: &str) -> bool {
        schema == BASE
    }

    fn type_definition(&self, schema: &str, kind: u16) -> Option<&TypeDefinition> {
        (schema == BASE)
            .then(|| self.definitions.get(&kind))
            .flatten()
    }

    fn lookup_type(&self, schema: &str, kind: u16) -> SchemaTypeLookup<'_> {
        if schema != BASE {
            return SchemaTypeLookup::Unknown;
        }
        if let Some(definition) = self.definitions.get(&kind) {
            return SchemaTypeLookup::Defined(definition);
        }
        if kind > BASE_TYPE_BOUND {
            SchemaTypeLookup::Absent
        } else {
            SchemaTypeLookup::PresentUnsupported
        }
    }

    fn supports_schema_key(&self, key: &SchemaKey) -> bool {
        key.raw() == ABC_KEY
    }
}

/// Types with explicitly established B-rep semantics; metadata retains its decoder name.
fn type_role(kind: u16) -> Option<&'static str> {
    Some(match kind {
        12 => "BODY",
        13 => "SHELL",
        14 => "FACE",
        15 => "LOOP",
        16 => "EDGE",
        17 => "HALFEDGE",
        18 => "VERTEX",
        19 => "REGION",
        29 => "POINT",
        30 => "LINE",
        31 => "CIRCLE",
        32 => "ELLIPSE",
        38 => "INTERSECTION",
        40 => "CHART",
        41 => "LIMIT",
        45 => "BSPLINE_VERTICES",
        50 => "PLANE",
        51 => "CYLINDER",
        52 => "CONE",
        53 => "SPHERE",
        124 => "B_SURFACE",
        126 => "NURBS_SURF",
        127 => "KNOT_MULT",
        128 => "KNOT_SET",
        133 => "TRIMMED_CURVE",
        134 => "B_CURVE",
        136 => "NURBS_CURVE",
        137 => "SP_CURVE",
        _ => return None,
    })
}

/// Adapter names are independent of input-declared extension names.
fn field_role(kind: u16, name: &str) -> Option<&'static str> {
    let geometry = matches!(kind, 30..=32 | 38 | 50..=53 | 124 | 133 | 134 | 137);
    if geometry {
        let common = match name {
            "local_id" => Some("node_id"),
            "owner_ref" => Some("owner"),
            "next_curve" | "next_surface" => Some("next"),
            "previous_curve" | "previous_surface" => Some("previous"),
            "indirect_owner" => Some("geometric_owner"),
            "orientation" => Some("sense"),
            _ => None,
        };
        if common.is_some() {
            return common;
        }
    }
    Some(match (kind, name) {
        (12, "size_precision") => "res_size",
        (12, "linear_precision") => "res_linear",
        (12, "body_kind") => "body_type",
        (12, "region_head") => "region",
        (12, "edge_head") => "edge",
        (12, "vertex_head") => "vertex",
        (13..=19 | 29, "local_id") => "node_id",
        (13, "next_shell")
        | (14, "next_back")
        | (15, "next_loop")
        | (16, "next_edge")
        | (18, "next_vertex")
        | (19, "next_region") => "next",
        (13, "back_faces") => "face",
        (13, "wire_edges") => "edge",
        (13, "isolated_vertex") => "vertex",
        (13, "region_ref") => "region",
        (13, "front_faces") => "front_face",
        (14, "loop_head") => "loop",
        (14, "back_shell") => "shell",
        (14, "surface_ref") => "surface",
        (14, "orientation") => "sense",
        (14, "front_shell") => "front_shell",
        (14, "next_front") => "next_front",
        (15 | 16, "fin_ref") => "halfedge",
        (15, "face_ref") => "face",
        (16 | 18, "precision") => "tolerance",
        (16, "curve_ref") => "curve",
        (16 | 18 | 29, "owner_ref") => "owner",
        (17, "loop_ref") => "loop",
        (17, "forward_fin") => "forward",
        (17, "backward_fin") => "backward",
        (17, "end_vertex") => "vertex",
        (17, "radial_fin") => "other",
        (17, "edge_ref") => "edge",
        (17, "curve_ref") => "curve",
        (17, "orientation") => "sense",
        (18, "point_ref") => "point",
        (19, "body_ref") => "body",
        (19, "shell_head") => "shell",
        (19, "region_kind") => "type",
        (29, "position") | (30 | 50..=52, "origin") => "pvec",
        (30, "tangent") => "direction",
        (31 | 32 | 53, "center") => "centre",
        (31 | 32 | 50..=53, "x_direction") => "x_axis",
        (31 | 32 | 50, "normal") => "normal",
        (31 | 51..=53, "radius") => "radius",
        (32, "major_radius") => "major_radius",
        (32, "minor_radius") => "minor_radius",
        (51..=53, "axis") => "axis",
        (52, "sin_half_angle") => "sin_half_angle",
        (52, "cos_half_angle") => "cos_half_angle",
        (38, "supporting_surfaces") => "surface",
        (38, "chart_ref") => "chart",
        (38, "start_limit") => "start",
        (38, "end_limit") => "end",
        (40, "base_parameter") => "base_parameter",
        (40, "base_scale") => "base_scale",
        (40, "chart_count") => "chart_count",
        (40, "chordal_error") => "chordal_error",
        (40, "angular_error") => "angular_error",
        (40, "sample_positions") | (41, "limit_positions") => "hvec",
        (41, "limit_kind") => "type",
        (124 | 134, "nurbs_ref") => "nurbs",
        (133, "basis_curve") => "basis_curve",
        (133, "start_point") => "point_1",
        (133, "end_point") => "point_2",
        (133, "start_parameter") => "parm_1",
        (133, "end_parameter") => "parm_2",
        (137, "surface_ref") => "surface",
        (137, "parameter_curve") => "b_curve",
        (137, "original_curve") => "original",
        (137, "original_tolerance") => "tolerance_to_original",
        (45, "vertices") => "vertices",
        (127, "multiplicities") => "mult",
        (128, "knots") => "knots",
        (136, "degree") => "degree",
        (136, "control_count") => "n_vertices",
        (126 | 136, "vertex_dimension") => "vertex_dim",
        (136, "knot_count") => "n_knots",
        (136, "knot_kind") => "knot_type",
        (136, "periodic") => "periodic",
        (136, "closed") => "closed",
        (126 | 136, "rational") => "rational",
        (136, "curve_form") => "curve_form",
        (126 | 136, "control_vertices") => "bspline_vertices",
        (136, "multiplicities") => "knot_mult",
        (136, "knots") => "knots",
        (126, "u_control_count") => "n_u_vertices",
        (126, "v_control_count") => "n_v_vertices",
        (126, "u_knot_count") => "n_u_knots",
        (126, "v_knot_count") => "n_v_knots",
        (126, "u_knot_kind") => "u_knot_type",
        (126, "v_knot_kind") => "v_knot_type",
        (126, "u_multiplicities") => "u_knot_mult",
        (126, "v_multiplicities") => "v_knot_mult",
        (126, "u_periodic") => "u_periodic",
        (126, "v_periodic") => "v_periodic",
        (126, "u_closed") => "u_closed",
        (126, "v_closed") => "v_closed",
        (126, "u_degree") => "u_degree",
        (126, "v_degree") => "v_degree",
        (126, "u_knots") => "u_knots",
        (126, "v_knots") => "v_knots",
        (126, "surface_form") => "surface_form",
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use parasolid_core::{FieldDefinition, FieldType};

    fn body_delta(insert_name: &str, replace_region: bool) -> SchemaResolution {
        let provider = AbcSchema::new().unwrap();
        let mut definition = provider.definitions[&12].clone();
        let original = definition.fields.clone();
        let extra = FieldDefinition {
            name: insert_name.to_owned(),
            field_type: FieldType::Integer,
            pointer_class: 0,
            element_count: 0,
            transmitted: true,
        };
        let mut edits = vec![SchemaEdit::Insert {
            offset: 0,
            field: extra.clone(),
        }];
        definition.fields = vec![extra];
        for field in original {
            if replace_region && field.name == "region" {
                edits.push(SchemaEdit::Delete { offset: 0 });
                edits.push(SchemaEdit::Insert {
                    offset: 0,
                    field: field.clone(),
                });
            } else {
                edits.push(SchemaEdit::Copy { offset: 0 });
            }
            definition.fields.push(field);
        }
        edits.push(SchemaEdit::End { offset: 0 });
        definition.source = SchemaSource::EmbeddedDelta;
        SchemaResolution {
            definition,
            edits,
            raw_schema: vec![],
            byte_range: 0..0,
        }
    }

    #[test]
    fn embedded_roles_require_copied_base_fields() {
        let provider = AbcSchema::new().unwrap();
        assert!(provider.audit(&[body_delta("metadata", false)]).is_ok());
        assert!(provider.audit(&[body_delta("region", false)]).is_err());
        assert!(provider.audit(&[body_delta("metadata", true)]).is_err());
    }

    #[test]
    fn partial_base_does_not_claim_unknown_types_are_absent() {
        let provider = AbcSchema::new().unwrap();
        assert!(matches!(
            provider.lookup_type(BASE, 12),
            SchemaTypeLookup::Defined(_)
        ));
        assert_eq!(
            provider.lookup_type(BASE, 100),
            SchemaTypeLookup::PresentUnsupported
        );
        assert_eq!(provider.lookup_type(BASE, 206), SchemaTypeLookup::Absent);
        assert_eq!(provider.lookup_type("13007", 12), SchemaTypeLookup::Unknown);
        assert!(provider.supports_schema_key(&SchemaKey::parse(ABC_KEY).unwrap()));
        assert!(
            !provider.supports_schema_key(&SchemaKey::parse("SCH_3001169_30100_13006").unwrap())
        );
    }
}
