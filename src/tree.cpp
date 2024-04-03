/*!
 * Inspired by the design of the tree in the xgboost and treelite package, both released under the Apache license 
 * with the following copyright:
 * Copyright 2015-2023 by XGBoost Contributors
 * Copyright 2017-2021 by [treelite] Contributors
 */
#include <stochtree/common.h>
#include <stochtree/tree.h>

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <sstream>

namespace StochTree {

constexpr std::int32_t Tree::kInvalidNodeId;
constexpr std::int32_t Tree::kDeletedNodeMarker;
constexpr std::int32_t Tree::kRoot;

std::int32_t Tree::NumLeaves() const {
  std::int32_t leaves { 0 };
  auto const& self = *this;
  this->WalkTree([&leaves, &self](std::int32_t nidx) {
                   if (self.IsLeaf(nidx)) {
                     leaves++;
                   }
                   return true;
                 });
  return leaves;
}

std::int32_t Tree::NumLeafParents() const {
  std::int32_t leaf_parents { 0 };
  auto const& self = *this;
  this->WalkTree([&leaf_parents, &self](std::int32_t nidx) {
                   if (!self.IsLeaf(nidx)){
                     if ((self.IsLeaf(self.LeftChild(nidx))) && (self.IsLeaf(self.RightChild(nidx)))){
                      leaf_parents++;
                     }
                   }
                   return true;
                 });
  return leaf_parents;
}

std::int32_t Tree::NumSplitNodes() const {
  std::int32_t splits { 0 };
  auto const& self = *this;
  this->WalkTree([&splits, &self](std::int32_t nidx) {
                   if (!self.IsLeaf(nidx)){
                     splits++;
                   }
                   return true;
                 });
  return splits;
}

void Tree::InplacePredictFromNodes(std::vector<double> result, std::vector<std::int32_t> node_indices) {
  if (result.size() != node_indices.size()) {
    Log::Fatal("Indices and result vector are different sizes");
  }
  data_size_t n = node_indices.size();
  for (data_size_t i = 0; i < n; i++) {
    result[i] = this->LeafValue(node_indices[i]);
  }
}

double Tree::PredictFromNode(std::int32_t node_id) {
  if (!this->IsLeaf(node_id)) {
    Log::Fatal("Node %d is not a leaf node", node_id);
  }
  return this->LeafValue(node_id);
}

std::vector<double> Tree::PredictFromNodes(std::vector<std::int32_t> node_indices) {
  data_size_t n = node_indices.size();
  std::vector<double> result(n);
  for (data_size_t i = 0; i < n; i++) {
    result[i] = PredictFromNode(node_indices[i]);
  }
  return result;
}

double Tree::PredictFromNode(std::int32_t node_id, Eigen::MatrixXd& basis, int row_idx) {
  if (!this->IsLeaf(node_id)) {
    Log::Fatal("Node %d is not a leaf node", node_id);
  }
  double pred = 0;
  for (int32_t k = 0; k < output_dimension_; k++) {
    pred += LeafValue(node_id, k) * basis(row_idx, k);
  }
  return pred;
}

std::vector<double> Tree::PredictFromNodes(std::vector<std::int32_t> node_indices, Eigen::MatrixXd& basis) {
  data_size_t n = node_indices.size();
  std::vector<double> result(n);
  for (data_size_t i = 0; i < n; i++) {
    result[i] = PredictFromNode(node_indices[i], basis, i);
  }
  return result;
}

void Tree::CloneFromTree(Tree* tree) {
  // Copy vectors from existing tree
  num_nodes = tree->num_nodes;
  num_deleted_nodes = tree->num_deleted_nodes;

  node_type_ = tree->node_type_;
  parent_ = tree->parent_;
  cleft_ = tree->cleft_;
  cright_ = tree->cright_;
  split_index_ = tree->split_index_;
  leaf_value_ = tree->leaf_value_;
  threshold_ = tree->threshold_;
  internal_nodes_ = tree->internal_nodes_;
  leaves_ = tree->leaves_;
  leaf_parents_ = tree->leaf_parents_;
  deleted_nodes_ = tree->deleted_nodes_;

  leaf_vector_ = tree->leaf_vector_;
  leaf_vector_begin_ = tree->leaf_vector_begin_;
  leaf_vector_end_ = tree->leaf_vector_end_;
  category_list_ = tree->category_list_;
  category_list_begin_ = tree->category_list_begin_;
  category_list_end_ = tree->category_list_end_;

  has_categorical_split_ = tree->has_categorical_split_;
  output_dimension_ = tree->output_dimension_;
}

std::int32_t Tree::AllocNode() {
  // Reuse a "deleted" node if available
  if (num_deleted_nodes != 0) {
    std::int32_t nid = deleted_nodes_.back();
    deleted_nodes_.pop_back();
    --num_deleted_nodes;
    return nid;
  }
  
  std::int32_t nd = num_nodes++;
  CHECK_LT(num_nodes, std::numeric_limits<int>::max());
  
  node_type_.push_back(TreeNodeType::kLeafNode);
  cleft_.push_back(kInvalidNodeId);
  cright_.push_back(kInvalidNodeId);
  split_index_.push_back(-1);
  leaf_value_.push_back(static_cast<double>(0));
  threshold_.push_back(static_cast<double>(0));
  // THIS is a placeholder, currently set after AllocNode is called ... 
  // ... to be refactored ...
  parent_.push_back(static_cast<double>(0));

  leaf_vector_begin_.push_back(leaf_vector_.size());
  leaf_vector_end_.push_back(leaf_vector_.size());
  category_list_begin_.push_back(category_list_.size());
  category_list_end_.push_back(category_list_.size());

  return nd;
}

void Tree::DeleteNode(std::int32_t nid) {
  CHECK_GE(nid, 1);
  auto pid = this->Parent(nid);
  bool is_left = this->LeftChild(pid) == nid;
  if (is_left) {
    SetLeftChild(pid, kInvalidNodeId);
  } else {
    SetRightChild(pid, kInvalidNodeId);
  }

  deleted_nodes_.push_back(nid);
  ++num_deleted_nodes;

  // Remove from vectors that track leaves, leaf parents, internal nodes, etc...
  leaves_.erase(std::remove(leaves_.begin(), leaves_.end(), nid), leaves_.end());
  leaf_parents_.erase(std::remove(leaf_parents_.begin(), leaf_parents_.end(), nid), leaf_parents_.end());
  internal_nodes_.erase(std::remove(internal_nodes_.begin(), internal_nodes_.end(), nid), internal_nodes_.end());
}

void Tree::ExpandNode(std::int32_t nid, int split_index, double split_value, double left_value, double right_value) {
  CHECK_EQ(output_dimension_, 1);
  int pleft = this->AllocNode();
  int pright = this->AllocNode();
  this->SetChildren(nid, pleft, pright);
  this->SetParents(nid, pleft, pright);
  this->SetNumericSplit(nid, split_index, split_value);
  this->SetLeaf(pleft, left_value);
  this->SetLeaf(pright, right_value);

  // Remove nid from leaves and add to internal nodes and leaf parents
  leaves_.erase(std::remove(leaves_.begin(), leaves_.end(), nid), leaves_.end());
  leaf_parents_.push_back(nid);
  internal_nodes_.push_back(nid);

  // Remove nid's parent node (if applicable) from leaf parents
  if (!IsRoot(nid)){
    std::int32_t parent_idx = Parent(nid);
    leaf_parents_.erase(std::remove(leaf_parents_.begin(), leaf_parents_.end(), parent_idx), leaf_parents_.end());
  }

  // Add pleft and pright to leaves
  leaves_.push_back(pleft);
  leaves_.push_back(pright);
}

void Tree::ExpandNode(std::int32_t nid, int split_index, std::vector<std::uint32_t> const& categorical_indices, double left_value, double right_value) {
  CHECK_EQ(output_dimension_, 1);
  int pleft = this->AllocNode();
  int pright = this->AllocNode();
  this->SetChildren(nid, pleft, pright);
  this->SetParents(nid, pleft, pright);
  this->SetCategoricalSplit(nid, split_index, categorical_indices);
  this->SetLeaf(pleft, left_value);
  this->SetLeaf(pright, right_value);

  // Remove nid from leaves and add to internal nodes and leaf parents
  leaves_.erase(std::remove(leaves_.begin(), leaves_.end(), nid), leaves_.end());
  leaf_parents_.push_back(nid);
  internal_nodes_.push_back(nid);

  // Remove nid's parent node (if applicable) from leaf parents
  if (!IsRoot(nid)){
    std::int32_t parent_idx = Parent(nid);
    leaf_parents_.erase(std::remove(leaf_parents_.begin(), leaf_parents_.end(), parent_idx), leaf_parents_.end());
  }

  // Add pleft and pright to leaves
  leaves_.push_back(pleft);
  leaves_.push_back(pright);
}

void Tree::ExpandNode(std::int32_t nid, int split_index, double split_value, std::vector<double> left_value_vector, std::vector<double> right_value_vector) {
  CHECK_GT(output_dimension_, 1);
  CHECK_EQ(output_dimension_, left_value_vector.size());
  CHECK_EQ(output_dimension_, right_value_vector.size());
  int pleft = this->AllocNode();
  int pright = this->AllocNode();
  this->SetChildren(nid, pleft, pright);
  this->SetParents(nid, pleft, pright);
  this->SetNumericSplit(nid, split_index, split_value);
  this->SetLeafVector(pleft, left_value_vector);
  this->SetLeafVector(pright, right_value_vector);

  // Remove nid from leaves and add to internal nodes and leaf parents
  leaves_.erase(std::remove(leaves_.begin(), leaves_.end(), nid), leaves_.end());
  leaf_parents_.push_back(nid);
  internal_nodes_.push_back(nid);

  // Remove nid's parent node (if applicable) from leaf parents
  if (!IsRoot(nid)){
    std::int32_t parent_idx = Parent(nid);
    leaf_parents_.erase(std::remove(leaf_parents_.begin(), leaf_parents_.end(), parent_idx), leaf_parents_.end());
  }

  // Add pleft and pright to leaves
  leaves_.push_back(pleft);
  leaves_.push_back(pright);
}

void Tree::ExpandNode(std::int32_t nid, int split_index, std::vector<std::uint32_t> const& categorical_indices, std::vector<double> left_value_vector, std::vector<double> right_value_vector) {
  CHECK_GT(output_dimension_, 1);
  CHECK_EQ(output_dimension_, left_value_vector.size());
  CHECK_EQ(output_dimension_, right_value_vector.size());
  int pleft = this->AllocNode();
  int pright = this->AllocNode();
  this->SetChildren(nid, pleft, pright);
  this->SetParents(nid, pleft, pright);
  this->SetCategoricalSplit(nid, split_index, categorical_indices);
  this->SetLeafVector(pleft, left_value_vector);
  this->SetLeafVector(pright, right_value_vector);

  // Remove nid from leaves and add to internal nodes and leaf parents
  leaves_.erase(std::remove(leaves_.begin(), leaves_.end(), nid), leaves_.end());
  leaf_parents_.push_back(nid);
  internal_nodes_.push_back(nid);

  // Remove nid's parent node (if applicable) from leaf parents
  if (!IsRoot(nid)){
    std::int32_t parent_idx = Parent(nid);
    leaf_parents_.erase(std::remove(leaf_parents_.begin(), leaf_parents_.end(), parent_idx), leaf_parents_.end());
  }

  // Add pleft and pright to leaves
  leaves_.push_back(pleft);
  leaves_.push_back(pright);
}

void Tree::ExpandNode(std::int32_t nid, int split_index, TreeSplit& split, double left_value, double right_value) {
  CHECK_EQ(output_dimension_, 1);
  if (split.NumericSplit()) {
    ExpandNode(nid, split_index, split.SplitValue(), left_value, right_value);
  } else {
    ExpandNode(nid, split_index, split.SplitCategories(), left_value, right_value);
  }
}

void Tree::ExpandNode(std::int32_t nid, int split_index, TreeSplit& split, std::vector<double> left_value_vector, std::vector<double> right_value_vector) {
  CHECK_GT(output_dimension_, 1);
  if (split.NumericSplit()) {
    ExpandNode(nid, split_index, split.SplitValue(), left_value_vector, right_value_vector);
  } else {
    ExpandNode(nid, split_index, split.SplitCategories(), left_value_vector, right_value_vector);
  }
}

void Tree::Reset() {
  // Clear all of the vectors that define the tree structure
  node_type_.clear();
  cleft_.clear();
  cright_.clear();
  split_index_.clear();
  leaf_value_.clear();
  threshold_.clear();
  parent_.clear();

  num_nodes = 0;
  has_categorical_split_ = false;

  leaf_vector_.clear();
  leaf_vector_begin_.clear();
  leaf_vector_end_.clear();
  category_list_.clear();
  category_list_begin_.clear();
  category_list_end_.clear();

  leaves_.clear();
  leaf_parents_.clear();
  internal_nodes_.clear();

  // Set bool / integer variables to default values
  num_nodes = 0;
  num_deleted_nodes = 0;
  has_categorical_split_ = false;
  output_dimension_ = 1;
}

void Tree::Init(std::int32_t output_dimension) {
  CHECK_GE(output_dimension, 1);

  // Clear all of the vectors that define the tree structure
  node_type_.clear();
  cleft_.clear();
  cright_.clear();
  split_index_.clear();
  leaf_value_.clear();
  threshold_.clear();
  parent_.clear();

  num_nodes = 0;
  has_categorical_split_ = false;

  leaf_vector_.clear();
  leaf_vector_begin_.clear();
  leaf_vector_end_.clear();
  category_list_.clear();
  category_list_begin_.clear();
  category_list_end_.clear();

  leaves_.clear();
  leaf_parents_.clear();
  internal_nodes_.clear();

  // Set output dimension
  output_dimension_ = output_dimension;

  // Allocate root node
  int rid = AllocNode();
  SetChildren(rid, kInvalidNodeId, kInvalidNodeId);
  SetParent(rid, kInvalidNodeId);
  if (output_dimension == 1) {
    this->SetLeaf(rid, 0.0);
  } else {
    this->SetLeafVector(rid, std::vector<double>(output_dimension, 0.));
  }

  // Add rid as a leaf node
  leaves_.push_back(rid);
}

void Tree::SetNumericSplit(std::int32_t nid, std::int32_t split_index, double threshold) {
  split_index_.at(nid) = split_index;
  threshold_.at(nid) = threshold;
  node_type_.at(nid) = TreeNodeType::kNumericalSplitNode;
}

void Tree::SetCategoricalSplit(std::int32_t nid, std::int32_t split_index, std::vector<std::uint32_t> const& category_list) {
  // CHECK(CategoryList(nid).empty());
  std::size_t const begin = category_list_.size();
  std::size_t const end = begin + category_list.size();
  category_list_.insert(category_list_.end(), category_list.begin(), category_list.end());
  category_list_begin_.at(nid) = begin;
  category_list_end_.at(nid) = end;

  split_index_.at(nid) = split_index;
  node_type_.at(nid) = TreeNodeType::kCategoricalSplitNode;

  has_categorical_split_ = true;
}

void Tree::SetLeaf(std::int32_t nid, double value) {
  CHECK_EQ(output_dimension_, 1);
  leaf_value_.at(nid) = value;
  cleft_.at(nid) = -1;
  cright_.at(nid) = -1;
  node_type_.at(nid) = TreeNodeType::kLeafNode;
}

void Tree::SetLeafVector(std::int32_t nid, std::vector<double> const& node_leaf_vector) {
  CHECK_GT(output_dimension_, 1);
  CHECK_EQ(output_dimension_, node_leaf_vector.size());
  if (HasLeafVector(nid)) {
    if (node_leaf_vector.size() != output_dimension_) {
      Log::Fatal("node_leaf_vector must be same size as the vector output dimension");
    }
    if (node_leaf_vector.size() != (leaf_vector_end_.at(nid) - leaf_vector_begin_.at(nid))) {
      Log::Fatal("Existing vector output is not the same size as node_leaf_vector");
    }
    std::size_t begin = leaf_vector_begin_.at(nid);
    std::size_t end = leaf_vector_end_.at(nid);
    std::size_t counter = 0;
    for (std::size_t i = begin; i < end; i++) {
      leaf_vector_[i] = node_leaf_vector[counter];
      counter++;
    }
  } else {
    std::size_t begin = leaf_vector_.size();
    std::size_t end = begin + node_leaf_vector.size();
    leaf_vector_.insert(leaf_vector_.end(), node_leaf_vector.begin(), node_leaf_vector.end());
    leaf_vector_begin_.at(nid) = begin;
    leaf_vector_end_.at(nid) = end;
  }

  split_index_.at(nid) = -1;
  cleft_.at(nid) = kInvalidNodeId;
  cright_.at(nid) = kInvalidNodeId;
  node_type_.at(nid) = TreeNodeType::kLeafNode;
}

void TreeNodeVectorsToJson(json& obj, Tree* tree) {
  // Initialize a map with names of the node vectors and empty json arrays
  std::map<std::string, json> tree_array_map;
  tree_array_map.emplace(std::pair("node_type", json::array()));
  tree_array_map.emplace(std::pair("parent", json::array()));
  tree_array_map.emplace(std::pair("left", json::array()));
  tree_array_map.emplace(std::pair("right", json::array()));
  tree_array_map.emplace(std::pair("split_index", json::array()));
  tree_array_map.emplace(std::pair("leaf_value", json::array()));
  tree_array_map.emplace(std::pair("threshold", json::array()));
  tree_array_map.emplace(std::pair("leaf_vector_begin", json::array()));
  tree_array_map.emplace(std::pair("leaf_vector_end", json::array()));
  tree_array_map.emplace(std::pair("category_list_begin", json::array()));
  tree_array_map.emplace(std::pair("category_list_end", json::array()));

  // Extract only the non-deleted nodes into tree_array_map
//  bool node_deleted;
  for (int i = 0; i < tree->NumNodes(); i++) {
//    node_deleted = (std::find(tree->deleted_nodes_.begin(), tree->deleted_nodes_.end(), i)
//                    != tree->deleted_nodes_.end());
//    if (!node_deleted) {
      tree_array_map["node_type"].emplace_back(static_cast<int>(tree->node_type_[i]));
      tree_array_map["parent"].emplace_back(tree->parent_[i]);
      tree_array_map["left"].emplace_back(tree->cleft_[i]);
      tree_array_map["right"].emplace_back(tree->cright_[i]);
      tree_array_map["split_index"].emplace_back(tree->split_index_[i]);
      tree_array_map["leaf_value"].emplace_back(tree->leaf_value_[i]);
      tree_array_map["threshold"].emplace_back(tree->threshold_[i]);
      tree_array_map["leaf_vector_begin"].emplace_back(static_cast<int>(tree->leaf_vector_begin_[i]));
      tree_array_map["leaf_vector_end"].emplace_back(static_cast<int>(tree->leaf_vector_end_[i]));
      tree_array_map["category_list_begin"].emplace_back(static_cast<int>(tree->category_list_begin_[i]));
      tree_array_map["category_list_end"].emplace_back(static_cast<int>(tree->category_list_end_[i]));
//    }
  }
  
  // Unpack the map into the reference JSON object
  for (auto& pair : tree_array_map) {
    obj.emplace(pair);
  }
}

void MultivariateLeafVectorToJson(json& obj, Tree* tree) {
  json vec = json::array();
  if (tree->leaf_vector_.size() > 0) {
//    bool node_deleted;
    for (int i = 0; i < tree->NumNodes(); i++) {
      // node_deleted = (std::find(tree->deleted_nodes_.begin(), tree->deleted_nodes_.end(), i)
      //                 != tree->deleted_nodes_.end());
      // if (!node_deleted) {
        if ((tree->leaf_vector_begin_[i] != 0) && (tree->leaf_vector_end_[i] != 0)) {
          for (uint64_t ind = tree->leaf_vector_begin_[i]; ind < tree->leaf_vector_end_[i]; ind++) {
            vec.emplace_back(tree->leaf_vector_[ind]);
          }
        }
      // }
    }
  }
  obj.emplace("leaf_vector", vec);
}

void SplitCategoryVectorToJson(json& obj, Tree* tree) {
  json vec = json::array();
  if (tree->leaf_vector_.size() > 0) {
//    bool node_deleted;
    for (int i = 0; i < tree->NumNodes(); i++) {
      // node_deleted = (std::find(tree->deleted_nodes_.begin(), tree->deleted_nodes_.end(), i)
      //                 != tree->deleted_nodes_.end());
      // if (!node_deleted) {
        if ((tree->category_list_begin_[i] != 0) && (tree->category_list_end_[i] != 0)) {
          for (uint64_t ind = tree->category_list_begin_[i]; ind < tree->category_list_end_[i]; ind++) {
            vec.emplace_back(static_cast<int>(tree->category_list_[ind]));
          }
        }
      // }
    }
  }
  obj.emplace("category_list", vec);
}

json Tree::to_json() {
  json result_obj;
  // Store the non-array fields in json
  result_obj.emplace("num_nodes", this->NumNodes());
  result_obj.emplace("num_deleted_nodes", this->NumDeletedNodes());
  result_obj.emplace("has_categorical_split", this->has_categorical_split_);
  result_obj.emplace("output_dimension", this->output_dimension_);

  // Unpack the array based fields
  TreeNodeVectorsToJson(result_obj, this);
  MultivariateLeafVectorToJson(result_obj, this);
  SplitCategoryVectorToJson(result_obj, this);
  
  // Initialize Json from Json::object map and return result
  return result_obj;
}

void JsonToTreeNodeVectors(const json& tree_json, Tree* tree) {
  int num_nodes = tree->NumNodes();
  for (int i = 0; i < num_nodes; i++) {
    tree->parent_.push_back(tree_json.at("parent").at(i));
    tree->cleft_.push_back(tree_json.at("left").at(i));
    tree->cright_.push_back(tree_json.at("right").at(i));
    tree->split_index_.push_back(tree_json.at("split_index").at(i));
    tree->leaf_value_.push_back(tree_json.at("leaf_value").at(i));
    tree->threshold_.push_back(tree_json.at("threshold").at(i));
    // Handle type conversions for node_type, leaf_vector_begin/end, and category_list_begin/end
    tree->node_type_.push_back(static_cast<TreeNodeType>(tree_json.at("node_type").at(i)));
    tree->leaf_vector_begin_.push_back(static_cast<uint64_t>(tree_json.at("leaf_vector_begin").at(i)));
    tree->leaf_vector_end_.push_back(static_cast<uint64_t>(tree_json.at("leaf_vector_end").at(i)));
    tree->category_list_begin_.push_back(static_cast<uint64_t>(tree_json.at("category_list_begin").at(i)));
    tree->category_list_end_.push_back(static_cast<uint64_t>(tree_json.at("category_list_end").at(i)));
  }
}

void JsonToMultivariateLeafVector(const json& tree_json, Tree* tree) {
  int num_entries = tree_json.at("leaf_vector").size();
  for (int i = 0; i < num_entries; i++) {
    tree->leaf_vector_.push_back(tree_json.at("leaf_vector").at(i));
  }
}

void JsonToSplitCategoryVector(const json& tree_json, Tree* tree) {
  int num_entries = tree_json.at("category_list").size();
  for (int i = 0; i < num_entries; i++) {
    tree->category_list_.push_back(tree_json.at("category_list").at(i));
  }
}

void Tree::from_json(const json& tree_json) {
  // Unpack non-array fields
  tree_json.at("num_nodes").get_to(this->num_nodes);
  tree_json.at("num_deleted_nodes").get_to(this->num_deleted_nodes);
  tree_json.at("has_categorical_split").get_to(this->has_categorical_split_);
  tree_json.at("output_dimension").get_to(this->output_dimension_);
  this->num_deleted_nodes = 0;
  
  // Unpack the array based fields
  JsonToTreeNodeVectors(tree_json, this);
  JsonToMultivariateLeafVector(tree_json, this);
  JsonToSplitCategoryVector(tree_json, this);
  
  // Reconstruct the internal_nodes, leaf_parents, and leaves vectors
  for (int i = 0; i < this->num_nodes; i++) {
    if (this->IsLeaf(i)) this->leaves_.push_back(i);
    else if (this->IsLeafParent(i)) this->leaf_parents_.push_back(i);
    else this->internal_nodes_.push_back(i);
  }
}

} // namespace StochTree
