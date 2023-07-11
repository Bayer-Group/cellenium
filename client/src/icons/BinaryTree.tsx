import Icon from "./svg/BinaryTree.svg";

interface Icon {
  size?: string;
  className?: string;
}

const BinaryTreeIcon = ({ className }: Icon) => {
  return <img src={Icon} className={className} alt="binary tree icon" />;
};

export default BinaryTreeIcon;
