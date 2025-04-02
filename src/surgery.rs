use crate::node::{FlatTree};



pub fn spr(
    flat_tree: &mut FlatTree,
    donor: usize,
    recipient: usize,
    time: f64,
) {
    let donor_parent = flat_tree[donor].parent.expect("The donor node should not be the root");
    let recipient_parent = flat_tree[recipient].parent.expect("The recipient node should not be the root");
    
    let recipient_sibling = if flat_tree[recipient_parent].left_child.unwrap() == recipient {
        flat_tree[recipient_parent].right_child.unwrap()
    } else {
        flat_tree[recipient_parent].left_child.unwrap()
    };

    let recipient_grandparent_opt = flat_tree[recipient_parent].parent;

    if donor_parent == recipient_parent {
        // In this case, only update the depth of the parent node of the donor (which is the same as the parent of the recipient).
        flat_tree[recipient_parent].depth = Some(time);
    } else {
        // Perform disconnects and reconnects to perform the SPR move. The nodes that must change are:
        // 1. The parent of the recipient node PR (changes father: it goes from PPR to PD, and changes child: it goes from SR to D).
        // 2. The sister of the recipient node SR (changes father: it goes from PD to PR).
        // 3. The parent of the donor node PD (changes child: it goes from D to PR).
        // 4. The donor node (changes parent: it goes from PD to PR).
        // 5. The grandparent of the recipient node PPR (changes child: it goes from PR to SR).

        // 1.
        // The root can change if the parent of the recipient is the root.
        // First, disconnect the receiver's parent, and give it the donor as a child. 
        // The receiver's parent has two children. When we disconnect PR, the sister of the receiver SR must be reconnected to the father of PR.
        // If there is 

        // If the recipient's parent is the root, the root of the tree changes. The recipient was one of the children of the root.
        // Now the root is the other child of the root.
        let recipient_has_grandparent =  !flat_tree[recipient_parent].parent.is_none();
        // The receiver's parent changes parent: the receiver's parent's parent is the donor's parent.
        flat_tree[recipient_parent].parent = Some(donor_parent);
        // The parent of the donor is now the child of the donor's parent.
        if flat_tree[recipient_parent].left_child.unwrap() == recipient {
            flat_tree[recipient_parent].left_child = Some(donor);
        } else {
            flat_tree[recipient_parent].right_child = Some(donor);
        }
        // Change the depth of the recipient's parent.
        flat_tree[recipient_parent].depth = Some(time);
        // Now the recipient's parent is correctly connected
        // It remains to modify: the receiver's parent's parent, the receiver's sister, the donor's parent, the donor
        // Fix the receiver's parent's parent: if the receiver was not the root, the receiver's parent's parent
        // now has a new child node, the sister of the receiver.
        // Also fix the sister of the receiver: it now has a new parent, the receiver's parent's parent.
        if recipient_has_grandparent {
            if flat_tree[recipient_grandparent_opt.unwrap()].left_child.unwrap() == recipient_parent {
                flat_tree[recipient_grandparent_opt.unwrap()].left_child = Some(recipient_sibling);
            } else {
                flat_tree[recipient_grandparent_opt.unwrap()].right_child = Some(recipient_sibling);
            }
            flat_tree[recipient_sibling].parent = recipient_grandparent_opt;
        } else {
            // If the recipient's parent was the root, the sister of the recipient becomes the new root.
            flat_tree[recipient_sibling].parent = None;
            flat_tree.root = recipient_sibling;
        }
        // Now the recipient's parent's parent and the recipient's sister are correctly connected.
        // We fix the donor's parent
        if flat_tree[donor_parent].left_child.unwrap() == donor {
            flat_tree[donor_parent].left_child = Some(recipient_parent);
        } else {
            flat_tree[donor_parent].right_child = Some(recipient_parent);
        }
        // The donor's parent is now correctly connected.
        // We fix the donor
        flat_tree[donor].parent = Some(recipient_parent);
        // The donor is now correctly connected.
        // The recipient's parent's parent, the recipient's sister, the donor's parent, and the donor are now correctly connected.
    }
}